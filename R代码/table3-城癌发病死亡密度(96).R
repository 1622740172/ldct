# ===================== 0. 环境彻底净化 =====================
# 【注释】标题，无实际运行效果，仅做代码标记
rm(list = ls())
# 【核心命令】删除R环境中的变量、数据、函数、对象
gc()
# 【辅助命令】垃圾回收，强制释放被删除对象占用的计算机内存

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
    
    # ===================== 【修改】随访时间转换（天转年） =====================
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
    as_lc_risk = as.numeric(as.character(as_lc_risk)),
    
    # 核心分组变量
    risk_group = case_when(
      as_lc_risk %in% c(1, 2) ~ "High",
      as_lc_risk == 0 ~ "Low",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(gender_label), !is.na(risk_group))

# 【修改】剔除发病和死亡时间均缺失的数据，保留至少有一项有效随访时间的记录
data_clean <- data_clean %>% filter(!is.na(cancer_years) | !is.na(death_years))

# ===================== 3. 核心计算函数 (泊松检验算 95% CI) =====================
# 【修改】增加 time_col 参数，根据不同事件动态选择时间列
calc_density_ci <- function(df, event_col, time_col) {
  # 仅针对该计算保留随访时间有效的人，防止有事件但无随访时间的情况
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

# 按照 Overall, Screened, Non-screened, Low 的顺序提取计算结果
generate_row_ci <- function(df, event_col, time_col) {
  res_overall <- calc_density_ci(df, event_col, time_col)
  res_screened <- calc_density_ci(df %>% filter(as_lc_risk == 1), event_col, time_col)
  res_nonscreened <- calc_density_ci(df %>% filter(as_lc_risk == 2), event_col, time_col)
  res_low <- calc_density_ci(df %>% filter(as_lc_risk == 0), event_col, time_col)
  
  return(c(res_overall["Events"], res_overall["PY"], res_overall["RateCI"],
           res_screened["Events"], res_screened["PY"], res_screened["RateCI"],
           res_nonscreened["Events"], res_nonscreened["PY"], res_nonscreened["RateCI"],
           res_low["Events"], res_low["PY"], res_low["RateCI"]))
}

# ===================== 4. 遍历所有协变量生成分层数据 =====================
cat("📈 正在按协变量进行分层计算人年及 95% CI...\n")

# 定义需要分层的协变量及其中文映射名（用于表格美化）
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

# 封装一个构建整个分层数据块的函数，增加 time_col
build_stratified_block <- function(data, event_col, time_col, title_name) {
  rows_list <- list()
  
  # 1. 插入大标题
  rows_list[[1]] <- c(title_name, rep(NA, 12))
  
  # 2. 插入总体 Total
  rows_list[[2]] <- c("Total", generate_row_ci(data, event_col, time_col))
  
  # 3. 遍历所有协变量
  for (cov in covariates) {
    if (cov %in% names(data)) {
      # 插入协变量标题
      rows_list[[length(rows_list) + 1]] <- c(cov_names[cov], rep(NA, 12))
      
      # 提取该协变量的各分类水平，排除NA
      levels_cov <- as.character(unique(data[[cov]]))
      levels_cov <- levels_cov[!is.na(levels_cov)]
      levels_cov <- sort(levels_cov) # 排个序更美观
      
      # 分层计算
      for (lvl in levels_cov) {
        subset_df <- data[data[[cov]] == lvl, ]
        row_data <- generate_row_ci(subset_df, event_col, time_col)
        # 用空格缩进表示层级关系
        rows_list[[length(rows_list) + 1]] <- c(paste0("  ", lvl), row_data)
      }
    }
  }
  
  # 转化为 DataFrame
  df_block <- do.call(rbind, rows_list)
  return(as.data.frame(df_block, stringsAsFactors = FALSE))
}

# ===================== 【修改】动态绑定对应的 time_col =====================
df_incidence <- build_stratified_block(data_clean, "LC_event", "cancer_years", "Lung cancer incidence")
df_mortality <- build_stratified_block(data_clean, "LC_death", "death_years", "Lung cancer mortality")
df_all_cause <- build_stratified_block(data_clean, "all_cause_death", "death_years", "All-cause mortality")

# 上下拼接，中间留空行分隔
empty_row <- as.data.frame(matrix(NA, nrow=1, ncol=13))
names(empty_row) <- names(df_incidence)
df_final <- rbind(df_incidence, empty_row, df_mortality, empty_row, df_all_cause)


# ===================== 5. Excel 复杂表头绘制与导出 =====================
out_dir <- "D:/HuaweiMoveData/Users/Ai/Desktop/非吸烟小论文"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_path <- file.path(out_dir, "全变量亚组发病密度分层表_最终版.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "Stratified Incidence")

# --- 画表头 ---
writeData(wb, "Stratified Incidence", "Overall", startRow = 1, startCol = 2)
writeData(wb, "Stratified Incidence", "High-risk group", startRow = 1, startCol = 5)
writeData(wb, "Stratified Incidence", "Low-risk group", startRow = 1, startCol = 11)
mergeCells(wb, "Stratified Incidence", cols = 2:4, rows = 1)
mergeCells(wb, "Stratified Incidence", cols = 5:10, rows = 1)
mergeCells(wb, "Stratified Incidence", cols = 11:13, rows = 1)

writeData(wb, "Stratified Incidence", "Screened", startRow = 2, startCol = 5)
writeData(wb, "Stratified Incidence", "Non-screened", startRow = 2, startCol = 8)
mergeCells(wb, "Stratified Incidence", cols = 5:7, rows = 2)
mergeCells(wb, "Stratified Incidence", cols = 8:10, rows = 2)

# 三级表头
sub_headers <- c("Events", "Person-years", "Rate (95% CI) per 100,000 person-years")
row3_headers <- c("Variables", rep(sub_headers, 4))
writeData(wb, "Stratified Incidence", t(row3_headers), startRow = 3, startCol = 1, colNames = FALSE)

# --- 写入拼接好的数据 ---
writeData(wb, "Stratified Incidence", df_final, startRow = 4, startCol = 1, colNames = FALSE)

# --- 格式美化 ---
center_style <- createStyle(halign = "center", valign = "center")
bold_style <- createStyle(textDecoration = "bold")
left_bold_style <- createStyle(halign = "left", valign = "center", textDecoration = "bold")
grey_bg_style <- createStyle(fgFill = "#F2F2F2", textDecoration = "bold")

addStyle(wb, "Stratified Incidence", center_style, rows = 1:3, cols = 1:13, gridExpand = TRUE)
addStyle(wb, "Stratified Incidence", bold_style, rows = 1:3, cols = 1:13, gridExpand = TRUE)

# ===================== 加粗3个大标题：发病、肺癌死亡、全因死亡 =====================
title_rows <- c(
  4,  # 肺癌发病
  4 + nrow(df_incidence) + 1,  # 肺癌死亡
  4 + nrow(df_incidence) + 1 + nrow(df_mortality) + 1  # 全因死亡
)
addStyle(wb, "Stratified Incidence", left_bold_style, rows = title_rows, cols = 1)

# 协变量行格式
cov_row_indices <- which(df_final[,1] %in% cov_names) + 3
if(length(cov_row_indices) > 0){
  for (r in cov_row_indices) {
    addStyle(wb, "Stratified Incidence", grey_bg_style, rows = r, cols = 1:13)
  }
}

# 调整列宽
setColWidths(wb, "Stratified Incidence", cols = 1, widths = 35)
setColWidths(wb, "Stratified Incidence", cols = c(2,5,8,11), widths = 10)
setColWidths(wb, "Stratified Incidence", cols = c(3,6,9,12), widths = 15)
setColWidths(wb, "Stratified Incidence", cols = c(4,7,10,13), widths = 38)

# 保存
tryCatch({
  saveWorkbook(wb, out_path, overwrite = TRUE)
  cat(paste0("\n✅ 发病/肺癌死亡/全因死亡 分层表生成成功！\n📁 文件已保存至：", out_path, "\n"))
}, error = function(e) {
  cat("\n❌ Excel保存失败：", e$message, "\n请确保文件没有被Excel打开占用。\n")
})