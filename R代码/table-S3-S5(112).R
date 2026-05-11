# ===================== 0. 彻底净化环境 =====================
# 清除工作区所有变量，防止被过往数据的残留污染
rm(list = ls())
# 回收系统内存
gc()

cat("\n✨ 环境已彻底净化，开始运行高危/低危基线特征对比代码...\n")

# ===================== 1. 加载包 =====================

library(dplyr)
library(tidyr)
library(data.table)
library(openxlsx)

# ===================== 2. 核心配置 =====================
file_path   <- "D:/匹配结果_最终修正版.csv"
# 输出路径文件名修改为高低危版
export_path <- "D:/HuaweiMoveData/Users/Ai/Desktop/112新版_肺癌患者高低危分布统计结果.xlsx"

# 筛查策略与分组变量设定
filter_cols <- c("Ncc_shai", "2024SHAI", "2024SHAI55-74")
filter_val  <- 1         
risk_col_name <- "as_lc_risk" # 【关键修改】高低危的分组列

age_breaks <- c(40, 45, 50, 55, 60, 65, 70, 75)
age_labels <- c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74")
bmi_breaks <- c(18.5, 24, 28)
bmi_labels <- c("<18.5", "18.5-24", "24-28", "≥28")

other_vars <- c(
  "as_occupational_hazard_exposure" = "职业暴露",
  "as_chronic_respiratory_disease"  = "慢性呼吸疾病", 
  "yijijiazushi"                    = "一级家族史",
  "all_second_smoke"                = "二手烟暴露"
)

# ===================== 3. 读取与清洗 =====================
cat("\n正在读取数据...\n")
df <- tryCatch({
  fread(file_path, encoding = "GBK", fill = TRUE, check.names = FALSE)
}, error = function(e) {
  fread(file_path, encoding = "UTF-8", fill = TRUE, check.names = FALSE)
})

# 处理重复列名
colnames(df) <- make.unique(colnames(df), sep = ".")
match_col <- function(target, all_cols) { grep(paste0("^", target, "(\\.|$)"), all_cols, value = TRUE)[1] }

age_col  <- match_col("age", colnames(df))
bmi_col  <- match_col("BMI", colnames(df))
risk_col <- match_col(risk_col_name, colnames(df)) # 匹配高低危列

filter_cols_matched <- sapply(filter_cols, match_col, colnames(df))
other_vars_matched  <- sapply(names(other_vars), match_col, colnames(df))

if (is.na(risk_col)) stop("⚠️ 找不到指定的高低危分组变量：", risk_col_name)

# 数据类型转换与高低危映射
df <- df %>%
  mutate(
    !!age_col := as.numeric(trimws(as.character(.data[[age_col]]))),
    !!bmi_col := as.numeric(trimws(as.character(.data[[bmi_col]]))),
    across(all_of(c(filter_cols_matched, other_vars_matched)), ~as.numeric(trimws(as.character(.x)))),
    # 【核心修改】将原来的 gender_clean 改为 risk_clean
    risk_clean = case_when(
      as.character(.data[[risk_col]]) %in% c("1", "2", "高危", "High", "high") ~ "高危",
      as.character(.data[[risk_col]]) %in% c("0", "低危", "Low", "low") ~ "低危",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(risk_clean)) # 过滤掉非高/低危的数据

# ===================== 4. 统计核心逻辑 =====================
# 标准化输出结构：筛选条件、变量名、变量分组、高危、低危、P值
standardize_output <- function(df, filter_name, var_name, group_name, p_val = NA) {
  if (!"高危" %in% names(df)) df$高危 <- "0 (0.0%)"
  if (!"低危" %in% names(df)) df$低危 <- "0 (0.0%)"
  
  df %>%
    mutate(
      筛选条件 = filter_name,
      变量名   = var_name,
      变量分组 = as.character(group_name),
      P        = as.character(p_val)
    ) %>%
    select(筛选条件, 变量名, 变量分组, 高危, 低危, P)
}

generate_example_table <- function(data, filter_col, filter_name, filter_val) {
  
  df_sub <- data %>% filter(.data[[filter_col]] == filter_val) %>% drop_na(all_of(c(age_col, "risk_clean")))
  if (nrow(df_sub) == 0) return(NULL)
  
  df_sub <- df_sub %>%
    mutate(
      age_group = factor(cut(.data[[age_col]], breaks = age_breaks, labels = age_labels, right = FALSE), levels = age_labels),
      BMI_group = factor(cut(.data[[bmi_col]], breaks = c(-Inf, bmi_breaks, Inf), labels = bmi_labels, right = FALSE), levels = bmi_labels)
    )
  
  # --- 1. 年龄均数 ---
  age_mean <- df_sub %>%
    group_by(risk_clean) %>%
    summarise(val = paste0(round(mean(.data[[age_col]], na.rm=T), 2), " ± ", round(sd(.data[[age_col]], na.rm=T), 2)), .groups = "drop") %>%
    pivot_wider(names_from = risk_clean, values_from = val) %>%
    standardize_output(filter_name, "年龄", "Mean ± SD")
  
  # --- 2. 年龄分组 ---
  age_group <- df_sub %>%
    count(age_group, risk_clean, .drop = FALSE) %>%
    group_by(risk_clean) %>%
    mutate(total = sum(n), val = ifelse(total==0, "0 (0.0%)", paste0(n, " (", round(n/total*100, 1), "%)"))) %>%
    ungroup() %>%
    pivot_wider(id_cols = age_group, names_from = risk_clean, values_from = val, values_fill = "0 (0.0%)")
  
  if (!"高危" %in% names(age_group)) age_group$高危 <- "0 (0.0%)"
  if (!"低危" %in% names(age_group)) age_group$低危 <- "0 (0.0%)"
  
  age_group <- age_group %>%
    mutate(筛选条件=filter_name, 变量名="年龄", 变量分组=as.character(age_group), P=NA_character_) %>%
    select(筛选条件, 变量名, 变量分组, 高危, 低危, P)
  
  # --- 3. BMI统计 ---
  bmi_stats <- df_sub %>%
    count(BMI_group, risk_clean, .drop = FALSE) %>%
    group_by(risk_clean) %>%
    mutate(total = sum(n), val = ifelse(total==0, "0 (0.0%)", paste0(n, " (", round(n/total*100, 1), "%)"))) %>%
    ungroup() %>%
    pivot_wider(id_cols = BMI_group, names_from = risk_clean, values_from = val, values_fill = "0 (0.0%)")
  
  if (!"高危" %in% names(bmi_stats)) bmi_stats$高危 <- "0 (0.0%)"
  if (!"低危" %in% names(bmi_stats)) bmi_stats$低危 <- "0 (0.0%)"
  
  bmi_tab <- table(df_sub$BMI_group, df_sub$risk_clean)
  bmi_p <- tryCatch({ if(any(chisq.test(bmi_tab)$expected < 5)) fisher.test(bmi_tab)$p.value else chisq.test(bmi_tab)$p.value }, error=function(e) NA)
  
  bmi_stats <- bmi_stats %>%
    mutate(筛选条件=filter_name, 变量名="BMI", 变量分组=as.character(BMI_group), 
           P = ifelse(row_number()==1, as.character(round(bmi_p, 3)), NA_character_)) %>%
    select(筛选条件, 变量名, 变量分组, 高危, 低危, P)
  
  # --- 4. 其他相关变量 ---
  other_res_list <- list()
  for (var_name in names(other_vars_matched)) {
    var_col <- other_vars_matched[var_name]
    var_label <- other_vars[var_name]
    
    temp_df <- df_sub %>%
      mutate(status = ifelse(.data[[var_col]] == 1, "是", "否")) %>%
      count(status, risk_clean) %>%
      group_by(risk_clean) %>%
      mutate(total = sum(n), val = ifelse(total==0, "0 (0.0%)", paste0(n, " (", round(n/total*100, 1), "%)"))) %>%
      ungroup() %>%
      pivot_wider(id_cols = status, names_from = risk_clean, values_from = val, values_fill = "0 (0.0%)")
    
    if (!"高危" %in% names(temp_df)) temp_df$高危 <- "0 (0.0%)"
    if (!"低危" %in% names(temp_df)) temp_df$低危 <- "0 (0.0%)"
    
    var_tab <- table(df_sub[[var_col]], df_sub$risk_clean)
    var_p <- tryCatch({ if(any(chisq.test(var_tab)$expected < 5)) fisher.test(var_tab)$p.value else chisq.test(var_tab)$p.value }, error=function(e) NA)
    
    temp_df <- temp_df %>%
      mutate(筛选条件=filter_name, 变量名=var_label, 变量分组=status,
             P = ifelse(row_number()==1, as.character(round(var_p, 3)), NA_character_)) %>%
      select(筛选条件, 变量名, 变量分组, 高危, 低危, P)
    
    other_res_list[[var_name]] <- temp_df
  }
  
  return(bind_rows(age_mean, age_group, bmi_stats, bind_rows(other_res_list)))
}

# ===================== 5. 执行循环与导出 =====================
cat("\n正在计算各策略基线特征...\n")
all_results <- list()
for (i in seq_along(filter_cols)) {
  f_col_m <- filter_cols_matched[i]
  if (is.na(f_col_m)) next
  res <- generate_example_table(df, f_col_m, filter_cols[i], filter_val)
  if (!is.null(res)) all_results[[filter_cols[i]]] <- res
}

final_df <- bind_rows(all_results)

# 最后检查并排序，确立标准列头
final_df <- final_df %>% select(筛选条件, 变量名, 变量分组, 高危, 低危, P)

# 输出至 Excel
wb <- createWorkbook()
addWorksheet(wb, "高低危基线特征")
header_style <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "TopBottom")
writeData(wb, "高低危基线特征", final_df, headerStyle = header_style)
setColWidths(wb, "高低危基线特征", cols = 1:ncol(final_df), widths = "auto")

saveWorkbook(wb, export_path, overwrite = TRUE)

cat("\n✅ 分析圆满完成！结果已成功导出至：\n", export_path, "\n")