# ===================== 0. 环境彻底净化 =====================
rm(
  list
  = ls())
gc()


# ===================== 0. 加载包 =====================
# 注意：已经不需要 survey 包了
library(data.table)
library(dplyr)
library(survival)
library(ggplot2)
library(tidyr)

# ===================== 1. 读取数据 =====================
file_path <- "D:/匹配结果_最终修正版.csv"
data <- fread(file_path)

# ===================== 2. 筛选与分组 =====================
# 选用原始策略 (as_lc_risk)
data_crude <- data %>%
  filter(as_lc_risk %in% c(1, 2)) %>%
  mutate(
    screening = ifelse(as_lc_risk == 1, 1, 0),
    # 【修改】：引入 cancer_time 和 death_time 取代 nian
    across(c(LC_event, LC_death, all_cause_death, cancer_time, death_time, age, BMI), as.numeric),
    # 【修改】：将天数转换为年
    cancer_years = cancer_time / 365.25,
    death_years = death_time / 365.25
  ) %>%
  mutate(
    # 保持之前定好的全部分组标准
    age_group = cut(age, breaks = c(-Inf, 50, 55, Inf), labels = c("40-49", "50-54", "55-74"), right = FALSE),
    BMI_group = cut(BMI, breaks = c(-Inf, 24, Inf), labels = c("<=24", ">24"), right = TRUE)
  )

# ===================== 3. 定义亚组变量与清洗 =====================
subgroup_vars <- c("gender", "age_group", "BMI_group", "all_second_smoke", 
                   "as_occupational_hazard_exposure", "as_chronic_respiratory_disease", 
                   "yijijiazushi")

data_crude <- data_crude %>% mutate(across(all_of(subgroup_vars), as.factor))

# 剔除缺失值
# 【修改】：从必填变量中剔除 nian
required_vars <- c("screening", subgroup_vars, "LC_event", "LC_death", "all_cause_death")
data_clean <- data_crude %>% 
  drop_na(all_of(required_vars)) %>% 
  # 【修改】：确保发病或死亡随访时间至少有一个有效且大于0
  filter((!is.na(cancer_years) & cancer_years > 0) | (!is.na(death_years) & death_years > 0))

cat("数据准备完成。当前分析样本量：", nrow(data_clean), "\n")
cat("【注意】：正在计算 Crude Rate Ratios，未使用任何加权或校正。\n")

# ===================== 4. 定义结局与循环框架 =====================
outcomes <- c("LC_event", "LC_death", "all_cause_death")
outcome_labels <- c("肺癌发病密度 (A)", "肺癌死亡率 (B)", "全因死亡率 (C)")

# 【新增】：建立结局指标与其对应的时间变量映射字典
time_vars <- c("LC_event" = "cancer_years",
               "LC_death" = "death_years",
               "all_cause_death" = "death_years")

covariate_levels <- lapply(subgroup_vars, function(v) levels(data_clean[[v]]))
names(covariate_levels) <- subgroup_vars

# ===================== 5. 单因素 (Crude) 亚组分析函数 =====================
# 【修改】：增加 time_var 参数
subgroup_analysis_crude <- function(data, outcome, time_var, covar, lvls) {
  res <- data.frame(Covariate = covar, Subgroup = lvls, 
                    HR = NA, Lower = NA, Upper = NA, P = NA, stringsAsFactors = FALSE)
  for (i in seq_along(lvls)) {
    sub_data <- data[data[[covar]] == lvls[i], ]
    
    # 【新增】：在子组内过滤掉当前结局随访时间无效的样本
    sub_data <- sub_data[!is.na(sub_data[[time_var]]) & sub_data[[time_var]] > 0, ]
    
    # 检查是否有两组数据且发生了事件
    if (nrow(sub_data) > 0 && length(unique(sub_data$screening)) == 2 && sum(sub_data[[outcome]]) > 0) {
      
      # 【核心改动】：最基础的单因素 Cox 模型，没有 weights，没有 design
      # 【修改】：动态使用对应的 time_var
      fit <- tryCatch({ 
        coxph(as.formula(paste0("Surv(", time_var, ", ", outcome, ") ~ screening")), data = sub_data) 
      }, error = function(e) NULL)
      
      if (!is.null(fit)) { 
        res[i, c("HR", "Lower", "Upper", "P")] <- c(
          exp(coef(fit))[1], 
          exp(confint(fit))[1, 1], 
          exp(confint(fit))[1, 2], 
          summary(fit)$coefficients[1, "Pr(>|z|)"]
        ) 
      }
    }
  }
  return(res)
}

# ===================== 6. 运行三大结局分析 =====================
all_results_crude <- list()

for (out_idx in seq_along(outcomes)) {
  out <- outcomes[out_idx]
  out_label <- outcome_labels[out_idx]
  current_time_var <- time_vars[out] # 【获取当前结局对应的时间变量】
  
  # 计算各亚组的 Crude HR
  res_df <- do.call(rbind, lapply(subgroup_vars, function(v) {
    # 【修改】：传入 current_time_var
    subgroup_analysis_crude(data_clean, out, current_time_var, v, covariate_levels[[v]])
  }))
  
  # 【修改】：计算总人群的 Crude HR 时，也要过滤当前事件对应的时间变量
  overall_data <- data_clean[!is.na(data_clean[[current_time_var]]) & data_clean[[current_time_var]] > 0, ]
  
  fit_all <- coxph(as.formula(paste0("Surv(", current_time_var, ", ", out, ") ~ screening")), data = overall_data)
  overall_row <- data.frame(
    Covariate = "Overall", 
    Subgroup = "Crude Total", 
    HR = exp(coef(fit_all))[1], 
    Lower = exp(confint(fit_all))[1, 1], 
    Upper = exp(confint(fit_all))[1, 2], 
    P = summary(fit_all)$coefficients[1, "Pr(>|z|)"], 
    stringsAsFactors = FALSE
  )
  
  all_results_crude[[out_label]] <- rbind(overall_row, res_df)
}

# ===================== 7. 绘制 Crude 森林图函数 =====================
plot_forest_crude <- function(df, out_name) {
  df$Covariate <- factor(df$Covariate, levels = c("Overall", subgroup_vars))
  df$Label <- ifelse(df$Covariate == "Overall", "Overall (Crude)", paste0(df$Covariate, ": ", df$Subgroup))
  
  df <- df %>% arrange(desc(Covariate == "Overall"), seq_len(nrow(df))) %>% 
    mutate(Label = factor(Label, levels = rev(Label)))
  
  # 生成数值标签
  df$Display_Text <- case_when(
    is.na(df$HR) | df$HR == 0 ~ "N/A", 
    TRUE ~ sprintf("%.2f (%.2f-%.2f)   P=%s", 
                   df$HR, df$Lower, df$Upper, 
                   ifelse(is.na(df$P), "-", ifelse(df$P < 0.001, "<0.001", sprintf("%.3f", df$P))))
  )
  
  ggplot(df, aes(x = HR, xmin = Lower, xmax = Upper, y = Label)) +
    geom_point(size = 2.5, color = "black") +  # 使用严谨的黑色表示未经校正的原始数据
    geom_errorbarh(height = 0.2, color = "black") + 
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") + # 1线标红更醒目
    geom_text(aes(x = 4, label = Display_Text), hjust = 0, size = 3.5, color = "black") + 
    scale_x_log10(name = "粗率比 / Crude Rate Ratio (95% CI)", breaks = c(0.1, 0.5, 1, 2)) + 
    coord_cartesian(xlim = c(0.1, 25)) + 
    labs(title = paste("原始策略 -", out_name), 
         subtitle = "未校正分析 (Unadjusted Subgroup Analysis)", 
         y = NULL) + 
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"), 
          panel.grid.minor = element_blank())
}

# ===================== 8. 输出图片与表格 =====================
# 【修改】：设置新的保存路径
output_dir <- "D:/HuaweiMoveData/Users/Ai/Desktop/非吸烟小论文/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (name in names(all_results_crude)) {
  # 提取干净的文件名后缀（去掉了标题里的A,B,C括号）
  file_suffix <- gsub(" \\(.*\\)", "", name) 
  
  p <- plot_forest_crude(all_results_crude[[name]], name)
  
  # 【修改】：增加路径前缀
  ggsave(paste0(output_dir, "原始策略_Crude森林图_", file_suffix, ".png"), p, width = 10.5, height = 7, dpi = 300)
  write.csv(all_results_crude[[name]], paste0(output_dir, "原始策略_Crude数据表_", file_suffix, ".csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

cat("\n✅ 所有 Crude Rate Ratios (粗率比) 分析运行完毕！图片与表格已保存至：", output_dir, "\n")