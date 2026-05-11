# ===================== 0. 环境彻底净化 =====================
rm(
  list
  = ls())
gc()


# ===================== 0. 加载包 =====================
library(data.table)
library(dplyr)
library(survival)
library(survey)
library(ggplot2)
library(tidyr)

# ===================== 1. 读取数据 =====================
file_path <- "D:/匹配结果_最终修正版.csv"
data <- fread(file_path)  # 自动检测编码

# ===================== 2. 筛选、预处理与强制分组 =====================
data <- data %>%
  filter(as_lc_risk %in% c(1, 2)) %>%
  mutate(
    screening = ifelse(as_lc_risk == 1, 1, 0),
    across(c(LC_event, LC_death, all_cause_death, cancer_time, death_time, age, BMI), as.numeric),
    cancer_years = cancer_time / 365.25,
    death_years = death_time / 365.25
  ) %>%
  mutate(
    # 【修改】：去掉 >=75 分组，如果存在 >=75 的人会被切成 NA 并在后续丢弃
    age_group = cut(age, 
                    breaks = c(-Inf, 50, 55, 75), 
                    labels = c("40-49", "50-54", "55-74"), 
                    right = FALSE),    
    # BMI分组 (right = TRUE 表示左开右闭，即 <=24 包含24本身)
    BMI_group = cut(BMI, 
                    breaks = c(-Inf, 24, Inf), 
                    labels = c("<=24", ">24"), 
                    right = TRUE)
  )

cat("数据清洗与分组完成。当前分析样本量：", nrow(data), "\n")

# ===================== 3. 定义协变量 (全部统一为分类变量) =====================
subgroup_vars <- c("gender", "age_group", "BMI_group", "all_second_smoke",
                   "as_occupational_hazard_exposure", "as_chronic_respiratory_disease",
                   "yijijiazushi")

# 将这些变量全部转为因子格式，确保后续回归不出错
data <- data %>% mutate(across(all_of(subgroup_vars), as.factor))

# ===================== 4. 剔除缺失值 =====================
required_vars <- c("screening", subgroup_vars, "LC_event", "LC_death", "all_cause_death")
data_clean <- data %>%
  drop_na(all_of(required_vars)) %>%
  filter((!is.na(cancer_years) & cancer_years > 0) | (!is.na(death_years) & death_years > 0))

cat("剔除缺失值后样本量：", nrow(data_clean), "\n")

# ===================== 5. 计算全局 IPW =====================
ps_formula <- as.formula(paste("screening ~", paste(subgroup_vars, collapse = " + ")))
ps_model <- glm(ps_formula, data = data_clean, family = binomial())
data_clean$ps <- predict(ps_model, type = "response")

data_clean$ipw <- ifelse(data_clean$screening == 1,
                         1 / data_clean$ps,
                         1 / (1 - data_clean$ps))

q <- quantile(data_clean$ipw, c(0.01, 0.99), na.rm = TRUE)
data_clean$ipw_trunc <- pmax(pmin(data_clean$ipw, q[2]), q[1])

# ===================== 6. 定义结局列表 =====================
outcomes <- c("LC_event", "LC_death", "all_cause_death")
outcome_labels <- c("肺癌发病密度", "肺癌死亡率", "全因死亡率")

time_vars <- c("LC_event" = "cancer_years",
               "LC_death" = "death_years",
               "all_cause_death" = "death_years")

covariate_levels <- list()
for (var in subgroup_vars) {
  covariate_levels[[var]] <- levels(data_clean[[var]])
}

# ===================== 7. 亚组分析函数 =====================
subgroup_analysis <- function(data, outcome, time_var, weight_var = "ipw_trunc",
                              covariate, levels, screening_var = "screening") {
  res <- data.frame(Covariate = covariate, Subgroup = levels, 
                    HR = NA, Lower = NA, Upper = NA, P = NA, stringsAsFactors = FALSE)
  for (i in seq_along(levels)) {
    sub_data <- data[data[[covariate]] == levels[i], ]
    
    sub_data <- sub_data[!is.na(sub_data[[time_var]]) & sub_data[[time_var]] > 0, ]
    
    if (nrow(sub_data) == 0) next
    if (length(unique(sub_data[[screening_var]])) < 2 || sum(sub_data[[outcome]]) == 0) next    
    
    design <- svydesign(ids = ~1, weights = ~get(weight_var), data = sub_data)
    formula <- as.formula(paste0("Surv(", time_var, ", ", outcome, ") ~ ", screening_var))
    
    fit <- tryCatch({ svycoxph(formula, design = design) }, error = function(e) NULL)    
    if (!is.null(fit)) {
      res[i, "HR"] <- exp(coef(fit))[1]
      res[i, "Lower"] <- exp(confint(fit))[1, 1]
      res[i, "Upper"] <- exp(confint(fit))[1, 2]
      res[i, "P"] <- summary(fit)$coefficients[1, "Pr(>|z|)"]
    }
  }
  return(res)
}

# ===================== 8. 执行分析 =====================
all_results <- list()

for (out_idx in seq_along(outcomes)) {
  outcome <- outcomes[out_idx]
  current_time_var <- time_vars[outcome]
  outcome_res <- data.frame()  
  
  for (var in subgroup_vars) {
    sub_res <- subgroup_analysis(data_clean, outcome, time_var = current_time_var, covariate = var, levels = covariate_levels[[var]])
    outcome_res <- rbind(outcome_res, sub_res)
  }  
  
  overall_data <- data_clean[!is.na(data_clean[[current_time_var]]) & data_clean[[current_time_var]] > 0, ]
  overall_design <- svydesign(ids = ~1, weights = ~ipw_trunc, data = overall_data)
  overall_formula <- as.formula(paste0("Surv(", current_time_var, ", ", outcome, ") ~ screening"))
  
  overall_fit <- svycoxph(overall_formula, design = overall_design)
  
  overall_row <- data.frame(
    Covariate = "Overall", Subgroup = "Total",
    HR = exp(coef(overall_fit))[1], Lower = exp(confint(overall_fit))[1, 1],
    Upper = exp(confint(overall_fit))[1, 2], P = summary(overall_fit)$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  outcome_res <- rbind(overall_row, outcome_res)
  all_results[[outcome_labels[out_idx]]] <- outcome_res
}

# ===================== 9. 绘图函数 =====================
plot_forest <- function(result_df, outcome_name) {  
  result_df$Covariate <- factor(result_df$Covariate, levels = c("Overall", subgroup_vars))
  result_df$Label <- ifelse(result_df$Covariate == "Overall", "Overall", paste0(result_df$Covariate, ": ", result_df$Subgroup))  
  result_df$Order <- seq_len(nrow(result_df))
  result_df <- result_df %>% arrange(desc(Covariate == "Overall"), Order)
  result_df$Label <- factor(result_df$Label, levels = rev(result_df$Label))  
  
  result_df <- result_df %>%
    mutate(
      Display_Text = case_when(
        is.na(HR) | HR == 0 ~ "N/A", 
        TRUE ~ sprintf("%.2f (%.2f-%.2f)   P=%s", 
                       HR, Lower, Upper, 
                       ifelse(is.na(P), "-", ifelse(P < 0.001, "<0.001", sprintf("%.3f", P))))
      )
    )  
  
  p <- ggplot(result_df, aes(x = HR, xmin = Lower, xmax = Upper, y = Label)) +
    geom_point(size = 2.5, color = "steelblue") +
    geom_errorbarh(height = 0.2, color = "steelblue") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    geom_text(aes(x = 4, label = Display_Text), hjust = 0, size = 3.5, color = "black", fontface = "plain") + 
    scale_x_log10(name = "校正后风险比 (95% CI)", breaks = c(0.1, 0.5, 1, 2)) +
    coord_cartesian(xlim = c(0.1, 25)) + 
    labs(title = paste("筛查有效性 -", outcome_name), subtitle = "亚组分析 (全分类变量)，IPW 加权", y = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 9),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          panel.grid.minor = element_blank())  
  return(p)
}

# ===================== 10. 结果输出与保存 =====================
# 【修改】：更新保存路径，并确保文件夹存在
output_dir <- "D:/HuaweiMoveData/Users/Ai/Desktop/非吸烟小论文/" 
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (outcome_name in names(all_results)) {
  p <- plot_forest(all_results[[outcome_name]], outcome_name)  
  
  # 保存图片
  ggsave(paste0(output_dir, "森林图_全分组版_", outcome_name, ".png"), 
         p, width = 10.5, height = 7, dpi = 300)  
  print(p)  
  
  # 保存CSV
  write.csv(all_results[[outcome_name]], 
            paste0(output_dir, "分类亚组分析_全分组版_", outcome_name, ".csv"), 
            row.names = FALSE, fileEncoding = "UTF-8")
}

cat("\n✅ 分析完成！文件已成功保存至：", output_dir, "\n")