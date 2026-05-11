# ===================== 0. 环境彻底净化 =====================
rm(list = ls())
gc()

# ===================== 1. 加载包 =====================
library(data.table)
library(dplyr)
library(survival)
library(survey)

# ===================== 2. 读取基础数据 =====================
file_path <- "D:/匹配结果_最终修正版.csv"
cat("🚀 正在读取数据...\n")
data <- tryCatch({
  fread(file_path, encoding = "GBK", quote = "", fill = TRUE, check.names = FALSE, data.table = FALSE)
}, error = function(e) {
  fread(file_path, encoding = "UTF-8", quote = "", fill = TRUE, check.names = FALSE, data.table = FALSE)
})



# ===================== 2. 读取与清洗基础数据 =====================
file_path <- "D:/匹配结果_最终修正版.csv"
cat("🚀 正在读取数据...\n")
data <- tryCatch({
  fread(file_path, encoding = "GBK", quote = "", fill = TRUE, check.names = FALSE, data.table = FALSE)
}, error = function(e) {
  fread(file_path, encoding = "UTF-8", quote = "", fill = TRUE, check.names = FALSE, data.table = FALSE)
})


# ---------------- 下方继续原来的数据清洗代码 ----------------
cat("正在清洗数据 (仅提取 as_lc_risk = 1 或 2 的原定高危人群)...\n")

data_clean <- data %>%
  filter(as_lc_risk %in% c(1, 2)) %>%
  mutate(
    screening = ifelse(as_lc_risk == 1, 1, 0),
    # ... 后面保持不变 ...

# ===================== 3. 数据清洗与 CanSPUC 高危人群提取 =====================

    
    # 格式化数值与时间
    BMI = as.numeric(as.character(BMI)),
    age = as.numeric(as.character(age)),
    cancer_time = as.numeric(as.character(cancer_time)),
    death_time = as.numeric(as.character(death_time)),
    cancer_years = cancer_time / 365.25, 
    death_years = death_time / 365.25,   
    
    LC_event = as.numeric(as.character(LC_event)),
    LC_death = as.numeric(as.character(LC_death)),
    all_cause_death = as.numeric(as.character(all_cause_death)),
    
    # 清洗性别与家族史变量，确保匹配准确
    gender = trimws(as.character(gender)),
    gender_label = case_when(
      gender == "M" ~ "男性",
      gender == "F" ~ "女性",
      TRUE ~ NA_character_
    ),
    yijijiazushi = as.character(yijijiazushi),
    
    # 生成基线协变量
    age_group = cut(age, breaks = c(-Inf, 50, 55, Inf), labels = c("40-49", "50-54", "55-74"), right = FALSE),
    BMI_group = cut(BMI, breaks = c(-Inf, 24, Inf), labels = c("<=24", ">24"), right = TRUE)
  )

# 定义需要用于 IPW 的 7 个协变量
covariates <- c("gender_label", "age_group", "BMI_group", "all_second_smoke", 
                "as_occupational_hazard_exposure", "as_chronic_respiratory_disease", 
                "yijijiazushi")

# 转换协变量为 factor 并剔除缺失值
data_clean <- data_clean %>%
  mutate(across(all_of(covariates), as.factor)) %>%
  tidyr::drop_na(all_of(c("screening", covariates, "LC_event", "LC_death", "all_cause_death")))

cat("✅ CanSPUC 高危人群清洗完毕，有效样本量：", nrow(data_clean), "\n")

# ===================== 4. 定义亚组分析核心函数 =====================
outcomes <- c("LC_event", "LC_death", "all_cause_death")
outcome_labels <- c("肺癌发病", "肺癌死亡", "全因死亡")
time_vars <- c("LC_event" = "cancer_years", 
               "LC_death" = "death_years", 
               "all_cause_death" = "death_years")

# 函数说明：在指定亚组内，剔除该分组变量重新计算 IPW，并跑三大结局 Cox
run_specific_subgroup <- function(data, sub_filter_expr, sub_label, exclude_var) {
  
  cat("\n--- 正在分析亚组:", sub_label, "---\n")
  
  # 1. 过滤出该亚组的数据
  sub_data <- data %>% filter(!!rlang::parse_expr(sub_filter_expr))
  
  if (nrow(sub_data) == 0) {
    cat("⚠️ 警告: 亚组数据为空，跳过分析。\n")
    return(NULL)
  }
  
  # 2. 重新计算该亚组内部的 IPW (剔除用来分组的变量，防止共线性报错)
  current_covs <- setdiff(covariates, exclude_var)
  ps_formula <- as.formula(paste("screening ~", paste(current_covs, collapse = " + ")))
  
  ps_model <- glm(ps_formula, data = sub_data, family = binomial())
  sub_data$ps <- predict(ps_model, type = "response")
  sub_data$ipw <- ifelse(sub_data$screening == 1, 1 / sub_data$ps, 1 / (1 - sub_data$ps))
  
  # 1% - 99% 截断
  q <- quantile(sub_data$ipw, c(0.01, 0.99), na.rm = TRUE)
  sub_data$ipw_trunc <- pmax(pmin(sub_data$ipw, q[2]), q[1])
  
  # 3. 运行三个主要结局的 Cox 回归
  results <- data.frame()
  
  for (out_idx in seq_along(outcomes)) {
    out <- outcomes[out_idx]
    out_lab <- outcome_labels[out_idx]
    t_var <- time_vars[out]
    
    # 过滤掉该结局随访时间 <= 0 的无效数据
    final_sub_data <- sub_data %>% filter(!is.na(.data[[t_var]]) & .data[[t_var]] > 0)
    
    if (nrow(final_sub_data) > 0 && length(unique(final_sub_data$screening)) == 2 && sum(final_sub_data[[out]]) > 0) {
      sub_design <- svydesign(ids = ~1, weights = ~ipw_trunc, data = final_sub_data)
      
      fit <- tryCatch({
        svycoxph(as.formula(paste0("Surv(", t_var, ", ", out, ") ~ screening")), design = sub_design)
      }, error = function(e) NULL)
      
      if (!is.null(fit)) {
        summ <- summary(fit)
        hr <- exp(coef(fit))[1]
        low <- exp(confint(fit))[1, 1]
        up <- exp(confint(fit))[1, 2]
        p_val <- summ$coefficients[1, "Pr(>|z|)"]
        
        res_row <- data.frame(
          Subgroup = sub_label,
          Outcome = out_lab,
          `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", hr, low, up),
          `P value` = ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val)),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        results <- rbind(results, res_row)
      }
    } else {
      # 如果事件数为0或只有单组
      res_row <- data.frame(
        Subgroup = sub_label, Outcome = out_lab, `HR (95% CI)` = "N/A", `P value` = "N/A", check.names = FALSE
      )
      results <- rbind(results, res_row)
    }
  }
  return(results)
}

# ===================== 5. 执行三大特定亚组分析 =====================

# ① 只保留女性
res_female <- run_specific_subgroup(data_clean, "gender_label == '女性'", "① 仅女性", "gender_label")

# ② 只保留年龄 55-74岁
res_age_5574 <- run_specific_subgroup(data_clean, "age_group == '55-74'", "② 年龄 55-74岁", "age_group")

# ③ 只保留有一级家族史
res_family <- run_specific_subgroup(data_clean, "yijijiazushi == '1'", "③ 有一级家族史", "yijijiazushi")

# ===================== 6. 结果汇总与导出 =====================
final_subgroup_results <- bind_rows(res_female, res_age_5574, res_family)

cat("\n================ 最终亚组分析结果 ================\n")
print(final_subgroup_results)

# 导出结果
out_dir <- "D:/HuaweiMoveData/Users/Ai/Desktop/非吸烟小论文/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_path <- paste0(out_dir, "CanSPUC_特定3大亚组主要结局分析.csv")

write.csv(final_subgroup_results, out_path, row.names = FALSE, fileEncoding = "GBK")
cat("\n🎉 分析完成！结果已保存至:", out_path, "\n")