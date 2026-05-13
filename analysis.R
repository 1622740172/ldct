# ============================================
# 亚组分析：比较筛查组与未筛查组的粗率比(RR)
# ============================================

# 安装必要的包
if (!require("data.table")) install.packages("data.table", repos="https://cloud.r-project.org/")
if (!require("ggplot2")) install.packages("ggplot2", repos="https://cloud.r-project.org/")
if (!require("gridExtra")) install.packages("gridExtra", repos="https://cloud.r-project.org/")

library(data.table)
library(ggplot2)
library(gridExtra)

# 读取数据
df <- fread("/workspace/匹配结果_最终修正版.csv", encoding="UTF-8")

# 筛选高危人群：as_lc_risk = 1 (筛查组) 或 2 (未筛查组)
df_high <- df[as_lc_risk %in% c(1, 2)]
df_high[, screened := ifelse(as_lc_risk == 1, 1, 0)]

# 计算人年
df_high[, cancer_time_years := cancer_time / 365.25]
df_high[, death_time_years := death_time / 365.25]

cat("高危人群总数:", nrow(df_high), "\n")
cat("筛查组:", sum(df_high$screened == 1), "\n")
cat("未筛查组:", sum(df_high$screened == 0), "\n")

# ============================================
# 计算粗率和RR的函数
# ============================================
calculate_rr <- function(dt, outcome_col, time_col) {
  screened <- dt[screened == 1]
  unscreened <- dt[screened == 0]
  
  e1 <- screened[, sum(get(outcome_col))]
  e0 <- unscreened[, sum(get(outcome_col))]
  py1 <- screened[, sum(get(time_col))]
  py0 <- unscreened[, sum(get(time_col))]
  
  rate1 <- ifelse(py1 > 0, (e1 / py1) * 1000, 0)
  rate0 <- ifelse(py0 > 0, (e0 / py0) * 1000, 0)
  rr <- ifelse(rate0 > 0, rate1 / rate0, NA)
  
  if (e1 > 0 && e0 > 0 && !is.na(rr) && !is.infinite(rr)) {
    log_rr <- log(rr)
    se_log_rr <- sqrt(1/e1 + 1/e0)
    ci_lower <- exp(log_rr - 1.96 * se_log_rr)
    ci_upper <- exp(log_rr + 1.96 * se_log_rr)
  } else {
    ci_lower <- NA
    ci_upper <- NA
  }
  
  return(list(
    screened_events = e1,
    unscreened_events = e0,
    screened_py = py1,
    unscreened_py = py0,
    screened_rate = rate1,
    unscreened_rate = rate0,
    rr = rr,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}

# ============================================
# 定义亚组结构
# ============================================
subgroups <- list(
  list(name = "Overall", var = NULL, val = NULL, is_header = FALSE),
  list(name = "Sex", var = "gender", val = "M", label = "Male", is_header = TRUE),
  list(name = "Sex", var = "gender", val = "F", label = "Female", is_header = FALSE),
  list(name = "Age group (years)", var = "age_group", val = 1, label = "40-49", is_header = TRUE),
  list(name = "Age group (years)", var = "age_group", val = 2, label = "50-54", is_header = FALSE),
  list(name = "Age group (years)", var = "age_group", val = 3, label = "55-74", is_header = FALSE),
  list(name = "Family history", var = "yijijiazushi", val = 0, label = "No", is_header = TRUE),
  list(name = "Family history", var = "yijijiazushi", val = 1, label = "Yes", is_header = FALSE),
  list(name = "Occupational hazard", var = "as_occupational_hazard_exposure", val = 0, label = "No", is_header = TRUE),
  list(name = "Occupational hazard", var = "as_occupational_hazard_exposure", val = 1, label = "Yes", is_header = FALSE),
  list(name = "Chronic respiratory disease", var = "as_chronic_respiratory_disease", val = 0, label = "No", is_header = TRUE),
  list(name = "Chronic respiratory disease", var = "as_chronic_respiratory_disease", val = 1, label = "Yes", is_header = FALSE)
)

outcomes <- list(
  LC_event = "cancer_time_years",
  LC_death = "death_time_years",
  all_cause_death = "death_time_years"
)

outcome_labels <- c(
  LC_event = "Lung cancer incidence",
  LC_death = "Lung cancer mortality", 
  all_cause_death = "All-cause mortality"
)

# ============================================
# 计算所有亚组的RR
# ============================================
all_results <- list()

for (outcome in names(outcomes)) {
  time_col <- outcomes[[outcome]]
  cat("\n========================================\n")
  cat("Outcome:", outcome, "\n")
  cat("========================================\n")
  
  res_list <- list()
  
  for (sg in subgroups) {
    if (is.null(sg$var)) {
      dt_sub <- df_high
      label <- "Overall"
    } else {
      dt_sub <- df_high[get(sg$var) == sg$val]
      label <- sg$label
    }
    
    res <- calculate_rr(dt_sub, outcome, time_col)
    
    res$outcome <- outcome
    res$group <- ifelse(is.null(sg$name), "Overall", sg$name)
    res$subgroup <- label
    res$is_header <- sg$is_header
    
    rr_str <- ifelse(is.na(res$rr) || is.infinite(res$rr), "N/A", sprintf("%.2f", res$rr))
    ci_str <- ifelse(is.na(res$ci_lower) || is.na(res$ci_upper), "N/A", 
                     sprintf("(%.2f-%.2f)", res$ci_lower, res$ci_upper))
    
    cat(sprintf("%s: RR=%s %s\n", label, rr_str, ci_str))
    
    res_list[[length(res_list) + 1]] <- res
  }
  
  all_results[[outcome]] <- res_list
}

# ============================================
# 准备绘图数据
# ============================================
plot_data <- data.frame()

for (outcome in names(outcomes)) {
  for (res in all_results[[outcome]]) {
    plot_data <- rbind(plot_data, data.frame(
      outcome = res$outcome,
      group = res$group,
      subgroup = res$subgroup,
      is_header = res$is_header,
      rr = res$rr,
      ci_lower = res$ci_lower,
      ci_upper = res$ci_upper,
      y_pos = 0
    ))
  }
}

# 为每个结局分配y位置
plot_data <- as.data.table(plot_data)
plot_data[, y_pos := .I, by = outcome]

# 反转y轴顺序（Overall在最上面）
max_y <- plot_data[, max(y_pos)]
plot_data[, y_pos := max_y - y_pos + 1]

# 标记分组标题行
plot_data[, is_group_header := FALSE]
plot_data[subgroup == "Overall", is_group_header := TRUE]
plot_data[is_header == TRUE & subgroup != "Overall", is_group_header := TRUE]

# 颜色定义
colors <- c(LC_event = "#2E86AB", LC_death = "#A23B72", all_cause_death = "#F18F01")

# ============================================
# 创建森林图函数
# ============================================
create_forest_plot <- function(dt, outcome_name, color) {
  dt_sub <- dt[outcome == outcome_name]
  
  # 创建标签
  dt_sub[, label := ifelse(is_group_header & subgroup != "Overall", 
                           paste0(group, "\n  ", subgroup),
                           ifelse(subgroup == "Overall", "Overall", paste0("  ", subgroup)))]
  
  # 创建RR文本
  dt_sub[, rr_text := ifelse(is.na(rr) | is.infinite(rr), "N/A", sprintf("%.2f", rr))]
  dt_sub[, ci_text := ifelse(is.na(ci_lower) | is.na(ci_upper), "N/A", 
                              sprintf("(%.2f-%.2f)", ci_lower, ci_upper))]
  dt_sub[, full_text := paste0(rr_text, "\n", ci_text)]
  
  # 过滤有效数据用于绘图
  dt_plot <- dt_sub[!is.na(rr) & !is.infinite(rr)]
  
  p <- ggplot() +
    # 置信区间
    geom_segment(data = dt_plot, 
                 aes(x = ci_lower, xend = ci_upper, y = y_pos, yend = y_pos),
                 color = color, linewidth = 1.2, alpha = 0.7) +
    # RR点
    geom_point(data = dt_plot,
               aes(x = rr, y = y_pos),
               color = color, size = 3.5, shape = 21, fill = color, stroke = 0.8) +
    # 参考线 RR=1
    geom_vline(xintercept = 1, color = "black", linewidth = 0.8) +
    # 对数刻度
    scale_x_log10(limits = c(0.15, 6),
                  breaks = c(0.5, 1, 1.5, 2, 2.5),
                  labels = c("0.5", "1", "1.5", "2", "2.5")) +
    scale_y_continuous(breaks = dt_sub$y_pos, labels = dt_sub$label) +
    coord_cartesian(ylim = c(0.5, max(dt_sub$y_pos) + 0.5)) +
    labs(title = outcome_labels[outcome_name], x = "RR (95% CI)", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(color = color, size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(10, 80, 10, 10)
    )
  
  # 添加数值标签
  p <- p + annotate("text", x = 4.5, y = dt_sub$y_pos, 
                    label = dt_sub$full_text, 
                    hjust = 0, vjust = 0.5, size = 3, color = color, fontface = "plain")
  
  return(p)
}

# 创建三个图
p1 <- create_forest_plot(plot_data, "LC_event", colors["LC_event"])
p2 <- create_forest_plot(plot_data, "LC_death", colors["LC_death"])
p3 <- create_forest_plot(plot_data, "all_cause_death", colors["all_cause_death"])

# 组合图
library(grid)
combined <- grid.arrange(
  p1, p2, p3,
  ncol = 3,
  top = textGrob(
    "Figure 1 The crude rate ratios of lung cancer incidence density, lung cancer mortality,\nand all-cause mortality comparing the screened with non-screened groups",
    gp = gpar(fontsize = 15, fontface = "bold")
  )
)

# 保存为TIFF格式 (300 dpi, 适合期刊投稿)
tiff("/workspace/forest_plot_R.tiff", width = 24, height = 14, units = "in", res = 300, compression = "lzw")
grid.draw(combined)
dev.off()

# 同时保存为PNG方便预览
png("/workspace/forest_plot_R.png", width = 2400, height = 1400, res = 150)
grid.draw(combined)
dev.off()

cat("\n图片已保存:\n")
cat("- /workspace/forest_plot_R.tiff (300 dpi, LZW压缩)\n")
cat("- /workspace/forest_plot_R.png (预览)\n")

# ============================================
# 导出结果表格
# ============================================
result_table <- data.frame()
for (outcome in names(outcomes)) {
  for (res in all_results[[outcome]]) {
    result_table <- rbind(result_table, data.frame(
      Outcome = outcome_labels[outcome],
      Group = res$group,
      Subgroup = res$subgroup,
      Screened_Events = res$screened_events,
      Screened_PY = round(res$screened_py, 1),
      Screened_Rate = round(res$screened_rate, 2),
      Unscreened_Events = res$unscreened_events,
      Unscreened_PY = round(res$unscreened_py, 1),
      Unscreened_Rate = round(res$unscreened_rate, 2),
      RR = ifelse(is.na(res$rr) || is.infinite(res$rr), NA, round(res$rr, 2)),
      CI_Lower = ifelse(is.na(res$ci_lower), NA, round(res$ci_lower, 2)),
      CI_Upper = ifelse(is.na(res$ci_upper), NA, round(res$ci_upper, 2))
    ))
  }
}

fwrite(result_table, "/workspace/rr_results_table.csv")
cat("\n结果表格已保存到 /workspace/rr_results_table.csv\n")
