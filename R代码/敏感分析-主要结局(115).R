# ===================== 0. 环境彻底净化 =====================
rm(list = ls())
gc()

# ===================== 1. 加载包 =====================
library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(survey)
library(openxlsx)

# ===================== 2. 读取与清洗基础数据 =====================
file_path <- "D:/匹配结果_最终修正版.csv"
cat("🚀 正在读取数据...\n")
data <- tryCatch({
  fread(file_path, encoding = "GBK", quote = "", fill = TRUE, check.names = FALSE, data.table = FALSE)
}, error = function(e) {
  fread(file_path, encoding = "UTF-8", quote = "", fill = TRUE, check.names = FALSE, data.table = FALSE)
})

cat("正在清洗数据 (仅提取 as_lc_risk = 1 或 2 的原定高危人群)...\n")

# 为了防止原始数据里未用到的大量重复列导致 dplyr 报错，
# 我们先单独把需要用到的核心列提取出来，再做清洗
core_cols <- c("as_lc_risk", "sensitive", "BMI", "age", "cancer_time", "death_time", 
               "LC_event", "LC_death", "all_cause_death", "gender", 
               "all_second_smoke", "as_occupational_hazard_exposure", 
               "as_chronic_respiratory_disease", "yijijiazushi")

# 提取核心列并清洗
data_clean <- data[, core_cols] %>%
  # 【核心回归】：不取交集！只按原始评估的高危/未高危进行筛选
  filter(as_lc_risk %in% c(1, 2)) %>%
  mutate(
    # 1为筛查组，2为未筛查组
    screening = ifelse(as_lc_risk == 1, 1, 0),
    
    # 敏感性分析排除标识
    sensitive = as.numeric(as.character(sensitive)),
    
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
    
    # 生成基线协变量
    age_group = cut(age, breaks = c(-Inf, 50, 55, Inf), labels = c("40-49", "50-54", "55-74"), right = FALSE),
    BMI_group = cut(BMI, breaks = c(-Inf, 24, Inf), labels = c("<=24", ">24"), right = TRUE),
    gender = as.factor(trimws(as.character(gender)))
  )

# 定义需要用于 IPW 的 7 个协变量
covariates <- c("gender", "age_group", "BMI_group", "all_second_smoke", 
                "as_occupational_hazard_exposure", "as_chronic_respiratory_disease", 
                "yijijiazushi")

# 转换协变量为 factor 并剔除缺失值
data_clean <- data_clean %>%
  mutate(across(all_of(covariates), as.factor)) %>%
  drop_na(all_of(c("screening", covariates, "LC_event", "LC_death", "all_cause_death")))

# ===================== 3. 定义 IPW 与 Cox 回归的核心函数 =====================
run_analysis <- function(df) {
  # 1. 重新计算该人群的 IPW 权重
  ps_formula <- as.formula(paste("screening ~", paste(covariates, collapse = " + ")))
  ps_model <- glm(ps_formula, data = df, family = binomial())
  df$ps <- predict(ps_model, type = "response")
  df$ipw <- ifelse(df$screening == 1, 1 / df$ps, 1 / (1 - df$ps))
  
  # 1% - 99% 截断
  q <- quantile(df$ipw, c(0.01, 0.99), na.rm = TRUE)
  df$ipw_trunc <- pmax(pmin(df$ipw, q[2]), q[1])
  
  # 2. 定义结局与时间映射
  outcomes <- c("LC_event", "LC_death", "all_cause_death")
  time_vars <- c("cancer_years", "death_years", "death_years")
  labels <- c("Lung cancer incidence", "Lung cancer mortality", "All-cause mortality")
  
  res_list <- list()
  
  # 3. 循环跑三个结局的加权 Cox
  for(i in 1:3) {
    out <- outcomes[i]
    t_var <- time_vars[i]
    
    # 过滤时间 <= 0 的无效数据
    sub_df <- df[!is.na(df[[t_var]]) & df[[t_var]] > 0, ]
    
    if(nrow(sub_df) > 0 && sum(sub_df[[out]]) > 0) {
      sub_design <- svydesign(ids = ~1, data = sub_df, weights = ~ipw_trunc)
      
      fit <- tryCatch({
        svycoxph(as.formula(paste0("Surv(", t_var, ", ", out, ") ~ screening")), design = sub_design)
      }, error = function(e) NULL)
      
      if(!is.null(fit)) {
        summ <- summary(fit)
        hr <- exp(coef(fit))[1]
        low <- exp(confint(fit))[1, 1]
        up <- exp(confint(fit))[1, 2]
        p_val <- summ$coefficients[1, "Pr(>|z|)"]
        
        res_list[[i]] <- data.frame(
          Outcomes = labels[i],
          HR_95CI = sprintf("%.2f (%.2f-%.2f)", hr, low, up),
          P_value = ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  return(do.call(rbind, res_list))
}

# ===================== 4. 分别运行主分析与敏感性分析 =====================
cat("📈 正在运行主分析 (包含全量 CanSPUC 高危人群)...\n")
res_primary <- run_analysis(data_clean)
colnames(res_primary)[2:3] <- c("Primary Analysis: HR (95% CI)", "Primary P value")

cat("📈 正在运行敏感性分析 (剔除半年内发病者)...\n")
# 敏感性分析数据：剔除 sensitive == 1 的人 
data_sensitive <- data_clean %>% filter(is.na(sensitive) | sensitive != 1)
res_sensitive <- run_analysis(data_sensitive)
colnames(res_sensitive)[2:3] <- c("Sensitivity Analysis*: HR (95% CI)", "Sensitivity P value")

# ===================== 5. 合并并美化表格 =====================
cat("🔗 正在合并对比表格...\n")
final_table <- left_join(res_primary, res_sensitive, by = "Outcomes")

# 重新排个序，让 HR 和 P 值挨着
final_table <- final_table %>%
  select(Outcomes, 
         `Primary Analysis: HR (95% CI)`, `Sensitivity Analysis*: HR (95% CI)`, 
         `Primary P value`, `Sensitivity P value`)

# ===================== 6. 导出到 Excel =====================
out_path <- "D:/HuaweiMoveData/Users/Ai/Desktop/非吸烟小论文/主效应敏感性分析对比表_CanSPUC人群_TableS6.xlsx"

wb <- createWorkbook()
addWorksheet(wb, "Sensitivity Analysis")

# 写数据
writeData(wb, "Sensitivity Analysis", final_table, startRow = 1)

# 写脚注
footnote <- "* Sensitivity analysis excluded participants diagnosed with lung cancer within 6 months of follow-up."
writeData(wb, "Sensitivity Analysis", footnote, startRow = nrow(final_table) + 3, startCol = 1)

# 格式美化
header_style <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "TopBottom")
addStyle(wb, "Sensitivity Analysis", header_style, rows = 1, cols = 1:ncol(final_table), gridExpand = TRUE)
setColWidths(wb, "Sensitivity Analysis", cols = 1, widths = 25)
setColWidths(wb, "Sensitivity Analysis", cols = 2:5, widths = 35)

saveWorkbook(wb, out_path, overwrite = TRUE)

cat("🎉 运算完成！敏感性分析对比表已生成。\n")
cat("📁 文件保存在：", out_path, "\n")