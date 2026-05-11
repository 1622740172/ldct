# ===================== 0. 环境彻底净化 =====================
rm(list = ls())
gc()

# ===================== 0. 加载必要包 =====================
library(data.table)
library(dplyr)
library(openxlsx)
library(stats) # 用于泊松精确检验

# ===================== 1. 读取数据 =====================
file_path <- "D:/匹配结果_最终修正版.csv"
cat("🚀 正在读取数据...\n")
data <- tryCatch({
  fread(file_path, encoding = "GBK", quote = "", fill = TRUE, check.names = FALSE, data.table = FALSE)
}, error = function(e) {
  fread(file_path, encoding = "UTF-8", quote = "", fill = TRUE, check.names = FALSE, data.table = FALSE)
})

# ===================== 2. 数据清洗与变量生成 =====================
cat("正在清洗数据及设定亚组变量...\n")
data_clean <- data %>%
  mutate(
    BMI = as.numeric(as.character(BMI)),
    
    # 随访时间转换（天转年）
    cancer_time = as.numeric(as.character(cancer_time)),
    death_time = as.numeric(as.character(death_time)),
    cancer_years = cancer_time / 365.25, # 肺癌发病随访时间 (年)
    death_years = death_time / 365.25,   # 死亡随访时间 (年)
    
    LC_event = as.numeric(as.character(LC_event)),
    LC_death = as.numeric(as.character(LC_death)),
    all_cause_death = as.numeric(as.character(all_cause_death)),
    
    BMI_GROUP = case_when(
      BMI < 24.0 ~ "<24.0",
      BMI >= 24.0 ~ "≥24.0",
      TRUE ~ NA_character_
    ),
    gender = trimws(as.character(gender)),
    gender_label = case_when(
      gender == "M" ~ "男性",
      gender == "F" ~ "女性",
      TRUE ~ NA_character_
    ),
    
    # 🌟 核心修改 1：提取策略列和新增的 as_lc_risk 变量
    as_lc_risk = as.numeric(as.character(.data[["as_lc_risk"]])),
    strategy_Ncc_shai = as.numeric(as.character(.data[["Ncc_shai"]])),
    strategy_2024SHAI = as.numeric(as.character(.data[["2024SHAI"]])),
    strategy_2024SHAI55_74 = as.numeric(as.character(.data[["2024SHAI55-74"]]))
  ) %>%
  filter(!is.na(gender_label)) %>%
  # 剔除发病和死亡时间均缺失的数据
  filter(!is.na(cancer_years) | !is.na(death_years))

# ===================== 3. 核心计算函数 (泊松检验算 95% CI) =====================
calc_density_ci <- function(df, event_col, time_col) {
  # 仅针对该计算保留随访时间有效的人
  valid_df <- df[!is.na(df[[time_col]]), ]
  
  events <- sum(valid_df[[event_col]], na.rm = TRUE)
  py_total <- sum(valid_df[[time_col]], na.rm = TRUE)
  
  if (py_total > 0) {
    ptest <- poisson.test(events, py_total)
    rate <- ptest$estimate * 100000
    ci_lower <- ptest$conf.int[1] * 100000
    ci_upper <- ptest$conf.int[2] * 100000
    rate_str <- sprintf("%.2f (%.2f-%.2f)", rate, ci_lower, ci_upper)
  } else {
    rate_str <- "-"
  }
  return(c(Events = as.character(events), PY = sprintf("%.2f", py_total), RateCI = rate_str))
}

# 🌟 核心修改 2：提取逻辑变为 Overall, High-screened, High-unscreened, Low
generate_row_ci <- function(df, event_col, time_col, strat_col) {
  res_overall <- calc_density_ci(df, event_col, time_col)
  
  # 高危筛查组 (策略=1 且 as_lc_risk=1)
  res_high_screened <- calc_density_ci(df %>% filter(!!sym(strat_col) == 1 & as_lc_risk == 1), event_col, time_col)
  
  # 高危未筛查组 (策略=1 且 as_lc_risk!=1，使用 is.na 防御空值报错)
  res_high_unscreened <- calc_density_ci(df %>% filter(!!sym(strat_col) == 1 & (is.na(as_lc_risk) | as_lc_risk != 1)), event_col, time_col)
  
  # 低危组 (策略=0)
  res_low <- calc_density_ci(df %>% filter(!!sym(strat_col) == 0), event_col, time_col)
  
  return(c(
    res_overall["Events"], res_overall["PY"], res_overall["RateCI"],
    res_high_screened["Events"], res_high_screened["PY"], res_high_screened["RateCI"],
    res_high_unscreened["Events"], res_high_unscreened["PY"], res_high_unscreened["RateCI"],
    res_low["Events"], res_low["PY"], res_low["RateCI"]
  ))
}

# ===================== 4. 定义生成整个分层表模块的函数 =====================
covariates <- c("gender_label", "age_group", "all_second_smoke", 
                "as_occupational_hazard_exposure", "as_chronic_respiratory_disease", 
                "yijijiazushi", "BMI_GROUP")

cov_names <- c("gender_label" = "Gender", 
               "age_group" = "Age group",
               "all_second_smoke" = "Second-hand smoke",
               "as_occupational_hazard_exposure" = "Occupational hazard exposure",
               "as_chronic_respiratory_disease" = "Chronic respiratory disease",
               "yijijiazushi" = "First-degree family history",
               "BMI_GROUP" = "BMI group")

build_stratified_block <- function(data, event_col, time_col, title_name, strat_col) {
  rows_list <- list()
  
  # 🌟 核心修改 3：大标题占位扩展至 13 列 (4组数据*3 + 1列表头)
  rows_list[[1]] <- c(title_name, rep(NA, 12))
  
  # 插入总体 Total
  rows_list[[2]] <- c("Total", generate_row_ci(data, event_col, time_col, strat_col))
  
  # 遍历所有协变量
  for (cov in covariates) {
    if (cov %in% names(data)) {
      rows_list[[length(rows_list) + 1]] <- c(cov_names[cov], rep(NA, 12))
      
      levels_cov <- as.character(unique(data[[cov]]))
      levels_cov <- levels_cov[!is.na(levels_cov)]
      levels_cov <- sort(levels_cov)
      
      for (lvl in levels_cov) {
        subset_df <- data[data[[cov]] == lvl, ]
        row_data <- generate_row_ci(subset_df, event_col, time_col, strat_col)
        rows_list[[length(rows_list) + 1]] <- c(paste0("  ", lvl), row_data)
      }
    }
  }
  
  df_block <- do.call(rbind, rows_list)
  return(as.data.frame(df_block, stringsAsFactors = FALSE))
}

# ===================== 5. 循环计算三种策略并存入列表 =====================
cat("📈 正在按策略计算发病/死亡密度及分层...\n")
strategies <- c("strategy_Ncc_shai", "strategy_2024SHAI", "strategy_2024SHAI55_74")
sheet_names <- c("Ncc_shai", "2024SHAI", "2024SHAI55-74")

results_list <- list()

for (i in seq_along(strategies)) {
  strat <- strategies[i]
  cat(paste0("  -> 正在处理策略: ", sheet_names[i], "\n"))
  
  df_inc <- build_stratified_block(data_clean, "LC_event", "cancer_years", "Lung cancer incidence", strat)
  df_mort <- build_stratified_block(data_clean, "LC_death", "death_years", "Lung cancer mortality", strat)
  df_all <- build_stratified_block(data_clean, "all_cause_death", "death_years", "All-cause mortality", strat)
  
  empty_row <- as.data.frame(matrix(NA, nrow=1, ncol=13)) # 扩展为13列
  names(empty_row) <- names(df_inc)
  
  df_final <- rbind(df_inc, empty_row, df_mort, empty_row, df_all)
  results_list[[sheet_names[i]]] <- df_final
}

# ===================== 6. Excel 复杂表头绘制与多Sheet导出 =====================
out_dir <- "D:/HuaweiMoveData/Users/Ai/Desktop/非吸烟小论文"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_path <- file.path(out_dir, "三种策略_全变量亚组发病死亡密度分层表.xlsx")

wb <- createWorkbook()

center_style <- createStyle(halign = "center", valign = "center")
bold_style <- createStyle(textDecoration = "bold")
left_bold_style <- createStyle(halign = "left", valign = "center", textDecoration = "bold")
grey_bg_style <- createStyle(fgFill = "#F2F2F2", textDecoration = "bold")

for (sheet_name in sheet_names) {
  addWorksheet(wb, sheet_name)
  df_target <- results_list[[sheet_name]]
  
  # 🌟 核心修改 4：四组一级表头坐标及合并更新
  writeData(wb, sheet_name, "Overall", startRow = 1, startCol = 2)
  writeData(wb, sheet_name, "High-risk screened", startRow = 1, startCol = 5)
  writeData(wb, sheet_name, "High-risk unscreened", startRow = 1, startCol = 8)
  writeData(wb, sheet_name, "Low-risk group", startRow = 1, startCol = 11)
  
  mergeCells(wb, sheet_name, cols = 2:4, rows = 1)
  mergeCells(wb, sheet_name, cols = 5:7, rows = 1)
  mergeCells(wb, sheet_name, cols = 8:10, rows = 1)
  mergeCells(wb, sheet_name, cols = 11:13, rows = 1)
  
  # 二级表头
  sub_headers <- c("Events", "Person-years", "Rate (95% CI) per 100,000 person-years")
  row2_headers <- c("Variables", rep(sub_headers, 4)) # 循环4次
  writeData(wb, sheet_name, t(row2_headers), startRow = 2, startCol = 1, colNames = FALSE)
  
  # 写入数据
  writeData(wb, sheet_name, df_target, startRow = 3, startCol = 1, colNames = FALSE)
  
  # 格式美化：扩展到 13 列
  addStyle(wb, sheet_name, center_style, rows = 1:2, cols = 1:13, gridExpand = TRUE)
  addStyle(wb, sheet_name, bold_style, rows = 1:2, cols = 1:13, gridExpand = TRUE)
  
  title_row_1 <- which(df_target[,1] == "Lung cancer incidence") + 2
  title_row_2 <- which(df_target[,1] == "Lung cancer mortality") + 2
  title_row_3 <- which(df_target[,1] == "All-cause mortality") + 2
  
  title_rows <- c(title_row_1, title_row_2, title_row_3)
  addStyle(wb, sheet_name, left_bold_style, rows = title_rows, cols = 1)
  
  cov_row_indices <- which(df_target[,1] %in% cov_names) + 2
  if(length(cov_row_indices) > 0){
    for (r in cov_row_indices) {
      addStyle(wb, sheet_name, grey_bg_style, rows = r, cols = 1:13)
    }
  }
  
  # 调整列宽 (适配 13 列)
  setColWidths(wb, sheet_name, cols = 1, widths = 35)
  setColWidths(wb, sheet_name, cols = c(2,5,8,11), widths = 10) # Events
  setColWidths(wb, sheet_name, cols = c(3,6,9,12), widths = 15) # Person-years
  setColWidths(wb, sheet_name, cols = c(4,7,10,13), widths = 38) # Rate
}

# 保存文件
tryCatch({
  saveWorkbook(wb, out_path, overwrite = TRUE)
  cat(paste0("\n✅ 密度计算与分层表生成完毕！\n📁 文件已保存至：", out_path, "\n"))
  cat("💡 提示：打开Excel文件，高危组已被成功分为 'screened' 和 'unscreened' 组。\n")
}, error = function(e) {
  cat("\n❌ Excel保存失败：", e$message, "\n请确保文件没有被Excel打开占用。\n")
})