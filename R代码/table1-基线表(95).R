# ===================== 0. 环境彻底净化 =====================
rm(list = ls())
gc()

# ===================== 0. 加载必要包 =====================
library(data.table)
library(dplyr)
library(tableone)
library(openxlsx)

# ===================== 1. 读取数据 =====================
file_path <- "D:/匹配结果_最终修正版.csv"
cat("🚀 正在读取数据...\n")
data <- tryCatch({
  fread(file_path, encoding = "GBK", quote = "", fill = TRUE, check.names = FALSE, data.table = FALSE)
}, error = function(e) {
  fread(file_path, encoding = "UTF-8", quote = "", fill = TRUE, check.names = FALSE, data.table = FALSE)
})

# ===================== 2. 数据清洗与变量生成 =====================
cat("正在生成 分组和风险分层变量...\n")
data_clean <- data %>%
  mutate(
    age = as.numeric(as.character(age)),
    BMI = as.numeric(as.character(BMI)),
    
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
    
    risk_group = case_when(
      as_lc_risk %in% c(1, 2) ~ "High",
      as_lc_risk == 0 ~ "Low",
      TRUE ~ NA_character_
    ),
    
    screen_group = case_when(
      as_lc_risk == 1 ~ "Screened",
      as_lc_risk == 2 ~ "Non-screened",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(gender_label), !is.na(risk_group))

# ===================== 3. 严格定义分析变量列表 =====================
vars_ordered <- c(
  "gender_label", 
  "age_group", "age",                                 
  "all_second_smoke", 
  "as_occupational_hazard_exposure", 
  "as_chronic_respiratory_disease",                   
  "yijijiazushi", 
  "BMI_GROUP", "BMI"                                  
)

factor_vars <- c("gender_label", "age_group", "all_second_smoke", 
                 "as_occupational_hazard_exposure", "as_chronic_respiratory_disease", 
                 "yijijiazushi", "BMI_GROUP")

existing_vars <- intersect(vars_ordered, names(data_clean))
existing_factors <- intersect(factor_vars, names(data_clean))

missing_vars <- setdiff(vars_ordered, names(data_clean))
if (length(missing_vars) > 0) {
  cat("⚠️ 警告：以下变量在数据中不存在：", paste(missing_vars, collapse = ", "), "\n")
}

data_clean <- data_clean %>% mutate(across(all_of(existing_factors), as.factor))

# ===================== 4. 生成双重指标并组合 =====================
cat("\n🚀 正在计算双重 SMD 和 P值，组装最终表头...\n")

# --- 计算 A：整体比对 (Total High vs Low) ---
tab_A <- CreateTableOne(vars = existing_vars, strata = "risk_group", data = data_clean, factorVars = existing_factors, addOverall = TRUE)
mat_A <- print(tab_A, smd = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE)
df_A <- as.data.frame(mat_A, stringsAsFactors = FALSE)

# --- 计算 B：高危内部比对 (Screened vs Non-screened) ---
data_high_only <- data_clean %>% filter(risk_group == "High")
tab_B <- CreateTableOne(vars = existing_vars, strata = "screen_group", data = data_high_only, factorVars = existing_factors, addOverall = FALSE)
mat_B <- print(tab_B, smd = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE)
df_B <- as.data.frame(mat_B, stringsAsFactors = FALSE)

# ===================== 5. 格式美化：Mean (SD) 转为 Mean ± SD =====================
format_continuous_row <- function(df, var_name) {
  row_idx <- which(rownames(df) == var_name)
  if (length(row_idx) > 0) {
    for (col in names(df)) {
      if (!col %in% c("SMD", "p", "test")) {
        df[row_idx, col] <- gsub(" \\(", " ± ", df[row_idx, col])
        df[row_idx, col] <- gsub("\\)", "", df[row_idx, col])
      }
    }
  }
  return(df)
}

if("age" %in% existing_vars) {
  df_A <- format_continuous_row(df_A, "age")
  df_B <- format_continuous_row(df_B, "age")
  rownames(df_A)[rownames(df_A) == "age"] <- "Age (mean ± SD)"
}
if("BMI" %in% existing_vars) {
  df_A <- format_continuous_row(df_A, "BMI")
  df_B <- format_continuous_row(df_B, "BMI")
  rownames(df_A)[rownames(df_A) == "BMI"] <- "BMI (mean ± SD)"
}

# --- 组装：加入P值列 ---
final_table <- data.frame(
  `Feature` = rownames(df_A),
  `overall` = df_A$Overall,
  `Total high-risk group` = df_A$High,
  `Screened` = df_B$Screened,
  `Non-screened` = df_B$`Non-screened`,
  `Standardized difference (screened vs non-screened)` = df_B$SMD,
  `P-value (screened vs non-screened)` = df_B$p,       # 🌟 新增：筛查与非筛查对比的P值
  `Low-risk group` = df_A$Low,
  `Standardized difference (high vs low)` = df_A$SMD,
  `P-value (high vs low)` = df_A$p,                    # 🌟 新增：高危与低危对比的P值
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(final_table) <- NULL

# ===================== 6. 稳健填充SMD和P值列 =====================
# 🌟 修改点：将P值列也加入到需要向下填充的列表中
fill_cols <- c(
  "Standardized difference (screened vs non-screened)", 
  "P-value (screened vs non-screened)",
  "Standardized difference (high vs low)",
  "P-value (high vs low)"
)

for(col in fill_cols) {
  final_table[[col]][trimws(final_table[[col]]) == ""] <- NA_character_
}

var_start_rows <- which(!grepl("^[[:space:]\\.]", final_table$Feature) & nchar(trimws(final_table$Feature)) > 0)

for (col in fill_cols) {
  for (i in seq_along(var_start_rows)) {
    start <- var_start_rows[i]
    end <- ifelse(i < length(var_start_rows), var_start_rows[i+1] - 1, nrow(final_table))
    
    fill_val <- NA_character_
    for (j in start:end) {
      if (!is.na(final_table[[col]][j])) {
        fill_val <- final_table[[col]][j]
        break
      }
    }
    if (!is.na(fill_val)) {
      final_table[[col]][start:end] <- fill_val
    }
  }
}

# ===================== 7. 导出Excel =====================
out_dir <- "D:/HuaweiMoveData/Users/Ai/Desktop/非吸烟小论文"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat(paste0("📁 已自动创建目标文件夹：", out_dir, "\n"))
}

out_path <- file.path(out_dir, "最终精简版_基线表_含P值与SMD.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "Baseline")
writeData(wb, "Baseline", final_table)

comment_text <- "说明：已精简变量，保留了慢呼史。分类变量格式为 n (%)，连续变量(Age, BMI)格式为 Mean ± SD。P值由卡方检验或T检验/方差分析得出。"
writeData(wb, "Baseline", comment_text, startRow = nrow(final_table) + 3, startCol = 1)

tryCatch({
  saveWorkbook(wb, out_path, overwrite = TRUE)
  cat(paste0("\n✅ 表格生成成功！\n📁 文件已保存至：", out_path, "\n"))
}, error = function(e) {
  cat("\n❌ Excel保存失败：", e$message, "\n")
  cat("建议：请确保文件没有被Excel打开占用。\n")
})