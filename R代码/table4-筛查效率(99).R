# ===================== 0. 环境彻底净化 =====================
rm(
  list
  = ls())
gc()


# ===================== 1. 环境准备 =====================

library(dplyr)
library(openxlsx)
library(data.table)

# ===================== 2. 核心配置 =====================
file_path     <- "D:/匹配结果_最终修正版.csv"
export_path   <- "D:/HuaweiMoveData/Users/Ai/Desktop/新数据/大论文结果/全人群-综合筛查效率指标表.xlsx"
strategy_cols <- c("as_lc_risk", "Ncc_shai", "2024SHAI", "2024SHAI55-74")
lc_event_col  <- "LC_event"
gender_col    <- "gender"
pass_values   <- list(
  as_lc_risk = c(1, 2),
  Ncc_shai = 1,
  `2024SHAI` = 1,
  `2024SHAI55-74` = 1
)
expected_rows <- 315008

# ===================== 3. 数据读取 =====================
cat("\n--- [1] 数据加载阶段 ---\n")
smart_read <- function(path, use_quote = "\"") {
  tryCatch({
    as.data.frame(fread(path, encoding = 'GBK', fill = TRUE, check.names = FALSE, quote = use_quote))
  }, error = function(e) {
    as.data.frame(fread(path, encoding = 'UTF-8', fill = TRUE, check.names = FALSE, quote = use_quote))
  })
}

df_full <- smart_read(file_path)
if (nrow(df_full) < expected_rows) df_full <- smart_read(file_path, use_quote = "")

actual_rows <- nrow(df_full)
cat("📊 实际读取行数:", actual_rows, "\n")

# ===================== 4. 数据清洗 =====================
# 统一转为数值型并处理空值
df_full[[lc_event_col]] <- as.numeric(trimws(as.character(df_full[[lc_event_col]])))
total_cancers <- sum(df_full[[lc_event_col]] == 1, na.rm = TRUE) # 全人群总肺癌数

# ===================== 5. 核心计算函数 (整合所有指标) =====================
calc_comprehensive_stats <- function(data, strat_col) {
  # A. 确定该策略的高危判定值
  pass_val <- if (strat_col %in% names(pass_values)) pass_values[[strat_col]] else 1
  
  # B. 标记高危人群 (High-risk)
  data_temp <- data %>%
    mutate(is_high_risk = ifelse(!!sym(strat_col) %in% pass_val, 1, 0))
  
  # C. 基础计数
  N_total     <- nrow(data_temp)                # 总样本量
  N_high_risk <- sum(data_temp$is_high_risk == 1, na.rm = TRUE) # 符合高危人数
  
  # D. 检出与漏诊计数
  # 检出：既是高危又是肺癌
  N_detected <- sum(data_temp$is_high_risk == 1 & data_temp[[lc_event_col]] == 1, na.rm = TRUE)
  # 漏诊：不是高危但确实是肺癌
  N_missed   <- sum(data_temp$is_high_risk == 0 & data_temp[[lc_event_col]] == 1, na.rm = TRUE)
  
  # E. 指标计算
  # 1. 总体符合率 (Overall Compliance/High-risk Rate)
  compliance_rate <- (N_high_risk / N_total) * 100
  
  # 2. 每10万人口符合LDCT筛查人数
  eligible_per_100k <- (N_high_risk / N_total) * 100000
  
  # 3. 每10万筛查者中检出的肺癌数
  detected_per_100k_screened <- if(N_high_risk > 0) (N_detected / N_high_risk) * 100000 else 0
  
  # 4. NNS (需筛查人数)
  nns_val <- if (N_detected > 0) N_high_risk / N_detected else Inf
  
  # 5. 漏诊率 (Missed case rate %)
  # 公式：基线非高危的新发肺癌 / 总新发肺癌
  missed_rate <- (N_missed / total_cancers) * 100
  
  return(tibble(
    筛查策略 = strat_col,
    符合高危人数 = N_high_risk,
    总体符合率_百分比 = paste0(round(compliance_rate, 2), "%"),
    每10万人口高危人数_Eligible = round(eligible_per_100k, 0),
    每10万筛查检出癌症数_Detected = round(detected_per_100k_screened, 1),
    需筛查人数_NNS = round(nns_val, 1),
    漏诊率_Missed_Rate_百分比 = paste0(round(missed_rate, 2), "%")
  ))
}

# ===================== 6. 执行计算与导出 =====================
cat("\n--- [2] 指标计算中 ---\n")

valid_strategy_cols <- intersect(strategy_cols, colnames(df_full))

final_results <- bind_rows(lapply(valid_strategy_cols, function(x) {
  calc_comprehensive_stats(df_full, x)
}))

# 打印预览
print(final_results)

# 写入 Excel
wb <- createWorkbook()
addWorksheet(wb, "综合筛查效率")
header_style <- createStyle(fontName = "Arial", textDecoration = "bold", 
                            fgFill = "#DCE6F1", border = "TopBottom", halign = "center")
writeData(wb, 1, final_results, headerStyle = header_style)
setColWidths(wb, 1, cols = 1:ncol(final_results), widths = 25)

# 保存
saveWorkbook(wb, export_path, overwrite = TRUE)

cat("\n✅ 处理完成！\n文件路径:", export_path, "\n")