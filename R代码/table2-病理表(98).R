# ===================== 0. 环境彻底净化 =====================
rm(list = ls())
gc()

# ===================== 0. 加载包 =====================
library(data.table)
library(dplyr)
library(openxlsx)

# ===================== 1. 文件路径和列名 =====================
file_path <- "D:/匹配结果_最终修正版.csv"   

stage_col   <- "fenqi"          # 病理分期
type_col    <- "bingli"         # 【修改】病理类型：0=不知道, 1=腺癌, 2=鳞癌, 3=小细胞癌, 4=其他
event_col   <- "LC_event"       # 1=肺癌发病
death_all   <- "all_cause_death"  # 全因死亡（所有人）
death_lc    <- "LC_death"       # 肺癌特异性死亡
risk_col    <- "as_lc_risk"     # 0=低风险, 1=筛查, 2=未筛查

# ===================== 2. 读取数据（解决GBK编码） =====================
data <- fread(file_path, encoding = "unknown", data.table = FALSE)

fix_encoding <- function(df, from = "GBK") {
  names(df) <- iconv(names(df), from = from, to = "UTF-8")
  for (col in names(df)) {
    if (is.character(df[[col]])) {
      df[[col]] <- iconv(df[[col]], from = from, to = "UTF-8", sub = "")
    }
  }
  return(df)
}
data <- fix_encoding(data)

# ===================== 3. 定义分组（总人群和肺癌患者） =====================
# 总人群（用于Total行、所有死亡统计的分母）
all_total <- data
hr_total  <- filter(data, get(risk_col) %in% c(1,2))
screened  <- filter(data, get(risk_col) == 1)
unscreened<- filter(data, get(risk_col) == 2)
low_risk  <- filter(data, get(risk_col) == 0)

# 肺癌患者（仅用于分期/类型统计）
all_lc   <- filter(data, get(event_col) == 1)
hr_lc    <- filter(all_lc, get(risk_col) %in% c(1,2))
scr_lc   <- filter(all_lc, get(risk_col) == 1)
unscr_lc <- filter(all_lc, get(risk_col) == 2)
low_lc   <- filter(all_lc, get(risk_col) == 0)

# 分母：总人群人数
N_total <- c(nrow(all_total), nrow(hr_total), nrow(screened), nrow(unscreened), nrow(low_risk))
# 分母：肺癌患者人数
N_lc    <- c(nrow(all_lc), nrow(hr_lc), nrow(scr_lc), nrow(unscr_lc), nrow(low_lc))

total_list <- list(all_total, hr_total, screened, unscreened, low_risk)
lc_list    <- list(all_lc,   hr_lc,   scr_lc,   unscr_lc,   low_lc)

# ===================== 4. 处理分期和组织学=====================
process_stage_histo <- function(df_lc) {
  df_lc %>%
    mutate(
      Stage = case_when(
        grepl("^I|^Ⅱ", get(stage_col), ignore.case = TRUE) ~ "I/II",
        grepl("^III|^Ⅳ", get(stage_col), ignore.case = TRUE) ~ "III/IV",
        TRUE ~ NA_character_
      ),
      # 【修改】严格匹配新的病理编码，空白或异常值归为 NA
      bingli_num = as.numeric(as.character(get(type_col))),
      Histo = case_when(
        bingli_num == 1 ~ "Adenocarcinoma",           # 1: 腺癌
        bingli_num == 2 ~ "Squamous cell carcinoma",  # 2: 鳞状细胞癌
        bingli_num == 3 ~ "Small cell carcinoma",     # 3: 小细胞癌
        bingli_num == 4 ~ "Others",                   # 4: 其他
        bingli_num == 0 ~ "Unknown",                  # 0: 不知道
        TRUE ~ NA_character_                          # 空白不用管
      )
    )
}

lc_list_processed <- lapply(lc_list, process_stage_histo)

# 已知分期的肺癌患者分母提取
lc_known_stage <- lapply(lc_list_processed, function(df) filter(df, !is.na(Stage)))
N_known_stage  <- sapply(lc_known_stage, nrow)

# 【新增】已知病理类型的肺癌患者分母提取（剔除了空白不用管的数据）
lc_known_histo <- lapply(lc_list_processed, function(df) filter(df, !is.na(Histo)))
N_known_histo  <- sapply(lc_known_histo, nrow)

# ===================== 5. 辅助函数 =====================
make_row <- function(df_list, n_list, var = NULL, value = NULL) {
  if (!is.null(var)) {
    counts <- sapply(df_list, function(df) sum(df[[var]] == value, na.rm = TRUE))
  } else {
    counts <- value
  }
  # 防止分母为0报错
  pct <- ifelse(n_list == 0, 0, round(counts / n_list * 100, 2))
  out <- paste0(counts, " (", pct, "%)")
  data.frame(t(out), stringsAsFactors = FALSE)
}

# ===================== 6. 逐行构建表格 =====================
col_names <- c("Feature", "overall", "Total high-risk group", 
               "Screened", "Non-screened", "Low-risk group")

rows <- list()

# ---- 1. Lung cancer incidence 标题 ----
rows[[length(rows)+1]] <- data.frame(Feature = "Lung cancer incidence", t(rep("-", 5)), stringsAsFactors = FALSE)
lc_counts <- sapply(total_list, function(df) sum(df[[event_col]] == 1, na.rm = TRUE))
total_row <- cbind(Feature = "Total", make_row(df_list = NULL, n_list = N_total, var = NULL, value = lc_counts))
rows[[length(rows)+1]] <- total_row

# ---- 2. Stage 标题 ----
rows[[length(rows)+1]] <- data.frame(Feature = "Stage", t(rep("-", 5)), stringsAsFactors = FALSE)
rows[[length(rows)+1]] <- cbind(Feature = "I/II", make_row(lc_known_stage, N_known_stage, var = "Stage", value = "I/II"))
rows[[length(rows)+1]] <- cbind(Feature = "III/IV", make_row(lc_known_stage, N_known_stage, var = "Stage", value = "III/IV"))

# ---- 3. Histological type 标题 【修改扩展分类】 ----
rows[[length(rows)+1]] <- data.frame(Feature = "Histological type", t(rep("-", 5)), stringsAsFactors = FALSE)
rows[[length(rows)+1]] <- cbind(Feature = "Adenocarcinoma",          make_row(lc_known_histo, N_known_histo, var = "Histo", value = "Adenocarcinoma"))
rows[[length(rows)+1]] <- cbind(Feature = "Squamous cell carcinoma", make_row(lc_known_histo, N_known_histo, var = "Histo", value = "Squamous cell carcinoma"))
rows[[length(rows)+1]] <- cbind(Feature = "Small cell carcinoma",    make_row(lc_known_histo, N_known_histo, var = "Histo", value = "Small cell carcinoma"))
rows[[length(rows)+1]] <- cbind(Feature = "Others",                  make_row(lc_known_histo, N_known_histo, var = "Histo", value = "Others"))
rows[[length(rows)+1]] <- cbind(Feature = "Unknown",                 make_row(lc_known_histo, N_known_histo, var = "Histo", value = "Unknown"))

# ---- 4. Deaths 标题 ----
rows[[length(rows)+1]] <- data.frame(Feature = "Deaths", t(rep("-", 5)), stringsAsFactors = FALSE)

# 全因死亡
all_cause_counts <- sapply(total_list, function(df) sum(df[[death_all]] == 1, na.rm = TRUE))
rows[[length(rows)+1]] <- cbind(Feature = "All-cause", make_row(df_list = NULL, n_list = N_total, var = NULL, value = all_cause_counts))

# 肺癌特异性死亡
lc_death_counts <- sapply(total_list, function(df) sum(df[[death_lc]] == 1, na.rm = TRUE))
rows[[length(rows)+1]] <- cbind(Feature = "Lung-cancer specific deaths", make_row(df_list = NULL, n_list = N_total, var = NULL, value = lc_death_counts))

# 合并所有行
final_table <- bind_rows(rows)
names(final_table) <- col_names

# ===================== 7. 导出 Excel =====================
out_path <- "D:/肺癌患者病理分布_最终修正版_细分病理.xlsx"
wb <- createWorkbook()
addWorksheet(wb, "Table")

writeData(wb, 1, final_table, startRow = 3)

# 双层表头
writeData(wb, 1, data.frame("", "overall", "High-risk group", "", "", "Low-risk group"), 
          startRow = 1, colNames = FALSE)
mergeCells(wb, 1, cols = 3:5, rows = 1)

writeData(wb, 1, data.frame("", "", "Total high-risk group", "Screened", "Non-screened", ""), 
          startRow = 2, colNames = FALSE)

setColWidths(wb, 1, cols = 1, widths = 35)
setColWidths(wb, 1, cols = 2:6, widths = 25)

# 样式美化
header_style <- createStyle(halign = "center", textDecoration = "bold", fontSize = 12)
addStyle(wb, 1, header_style, rows = 1:2, cols = 1:6, gridExpand = TRUE)

title_rows <- which(final_table$Feature %in% c("Lung cancer incidence", "Stage", 
                                               "Histological type", "Deaths")) + 2
title_style <- createStyle(textDecoration = "bold", halign = "center")
for (r in title_rows) addStyle(wb, 1, title_style, rows = r, cols = 1:6, gridExpand = TRUE)

data_style <- createStyle(halign = "center")
addStyle(wb, 1, data_style, rows = 3:(nrow(final_table)+2), cols = 1:6, gridExpand = TRUE)

saveWorkbook(wb, out_path, overwrite = TRUE)
cat("\n✅ 病理分期/分型表格已生成！已应用新的病理分型规则。\n")
cat("📁 文件保存在：", out_path, "\n")