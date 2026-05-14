# ============================================
# Figure 2: Rate Difference (RD) 柱状图
# 每个亚组显示三个结局的RD值（带误差线）
# ============================================

library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)

# 读取数据
df <- fread("/workspace/匹配结果_最终修正版.csv", encoding="UTF-8")

# 筛选高危人群
df_high <- df[as_lc_risk %in% c(1, 2)]
df_high[, screened := ifelse(as_lc_risk == 1, 1, 0)]

# 计算人年
df_high[, cancer_time_years := cancer_time / 365.25]
df_high[, death_time_years := death_time / 365.25]

# ============================================
# 计算Rate Difference
# ============================================
calculate_rd <- function(dt, outcome_col, time_col) {
  screened <- dt[screened == 1]
  unscreened <- dt[screened == 0]
  
  n1 <- nrow(screened)
  n0 <- nrow(unscreened)
  e1 <- screened[, sum(get(outcome_col))]
  e0 <- unscreened[, sum(get(outcome_col))]
  py1 <- screened[, sum(get(time_col))]
  py0 <- unscreened[, sum(get(time_col))]
  
  rate1 <- ifelse(py1 > 0, (e1 / py1) * 1000, 0)
  rate0 <- ifelse(py0 > 0, (e0 / py0) * 1000, 0)
  rd <- rate1 - rate0
  
  if (py1 > 0 && py0 > 0) {
    se_rd <- sqrt(e1 / py1^2 + e0 / py0^2) * 1000
    ci_lower <- rd - 1.96 * se_rd
    ci_upper <- rd + 1.96 * se_rd
  } else {
    se_rd <- NA
    ci_lower <- NA
    ci_upper <- NA
  }
  
  return(list(rd = rd, ci_lower = ci_lower, ci_upper = ci_upper))
}

# ============================================
# 定义亚组结构
# ============================================
subgroups <- list(
  list(type = "data", name = "Overall", var = NULL, val = NULL, label = "Overall"),
  list(type = "header", name = "Sex"),
  list(type = "data", name = "Sex", var = "gender", val = "M", label = "Male"),
  list(type = "data", name = "Sex", var = "gender", val = "F", label = "Female"),
  list(type = "header", name = "Age group (years)"),
  list(type = "data", name = "Age group (years)", var = "age_group", val = 1, label = "40-49"),
  list(type = "data", name = "Age group (years)", var = "age_group", val = 2, label = "50-54"),
  list(type = "data", name = "Age group (years)", var = "age_group", val = 3, label = "55-74"),
  list(type = "header", name = "Family history"),
  list(type = "data", name = "Family history", var = "yijijiazushi", val = 0, label = "No"),
  list(type = "data", name = "Family history", var = "yijijiazushi", val = 1, label = "Yes"),
  list(type = "header", name = "Occupational hazard"),
  list(type = "data", name = "Occupational hazard", var = "as_occupational_hazard_exposure", val = 0, label = "No"),
  list(type = "data", name = "Occupational hazard", var = "as_occupational_hazard_exposure", val = 1, label = "Yes"),
  list(type = "header", name = "Chronic respiratory disease"),
  list(type = "data", name = "Chronic respiratory disease", var = "as_chronic_respiratory_disease", val = 0, label = "No"),
  list(type = "data", name = "Chronic respiratory disease", var = "as_chronic_respiratory_disease", val = 1, label = "Yes")
)

outcomes <- list(
  LC_event = "cancer_time_years",
  LC_death = "death_time_years",
  all_cause_death = "death_time_years"
)

outcome_labels <- c(
  LC_event = "Lung cancer\nincidence",
  LC_death = "Lung cancer\nmortality", 
  all_cause_death = "All-cause\nmortality"
)

outcome_colors <- c(LC_event = "#2E86AB", LC_death = "#A23B72", all_cause_death = "#F18F01")

# ============================================
# 计算所有RD结果并整理为长格式
# ============================================
bar_data <- data.frame()

for (outcome in names(outcomes)) {
  time_col <- outcomes[[outcome]]
  
  for (sg in subgroups) {
    if (sg$type == "header") next
    
    if (is.null(sg$var)) {
      dt_sub <- df_high
      label <- "Overall"
      group_name <- "Overall"
    } else {
      dt_sub <- df_high[get(sg$var) == sg$val]
      label <- sg$label
      group_name <- sg$name
    }
    
    res <- calculate_rd(dt_sub, outcome, time_col)
    
    bar_data <- rbind(bar_data, data.frame(
      subgroup = label,
      group = group_name,
      outcome = outcome,
      outcome_label = outcome_labels[outcome],
      rd = res$rd,
      ci_lower = res$ci_lower,
      ci_upper = res$ci_upper,
      stringsAsFactors = FALSE
    ))
  }
}

# 设置亚组顺序
subgroup_order <- c("Overall", "Male", "Female", "40-49", "50-54", "55-74",
                    "No", "Yes", "No", "Yes", "No", "Yes")
group_order <- c("Overall", "Sex", "Sex", "Age group (years)", "Age group (years)", "Age group (years)",
                 "Family history", "Family history", "Occupational hazard", "Occupational hazard",
                 "Chronic respiratory disease", "Chronic respiratory disease")

# 创建唯一的亚组标识
bar_data$subgroup_id <- paste0(bar_data$group, "_", bar_data$subgroup)
unique_subgroups <- unique(bar_data$subgroup_id)

# 设置因子顺序
bar_data$subgroup_id <- factor(bar_data$subgroup_id, levels = unique_subgroups)
bar_data$outcome <- factor(bar_data$outcome, levels = names(outcomes))

# ============================================
# 创建分组柱状图
# ============================================

# 计算Y轴范围
y_min <- min(bar_data$ci_lower, na.rm = TRUE) * 1.1
y_max <- max(bar_data$ci_upper, na.rm = TRUE) * 1.1

# 创建分组标签
group_labels <- data.frame()
for (g in unique(bar_data$group)) {
  group_labels <- rbind(group_labels, data.frame(
    group = g,
    x = mean(which(bar_data$group[bar_data$outcome == "LC_event"] == g)),
    stringsAsFactors = FALSE
  ))
}

# 绘制柱状图
p <- ggplot(bar_data, aes(x = subgroup_id, y = rd, fill = outcome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                position = position_dodge(width = 0.8), width = 0.3, linewidth = 0.6, color = "black") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  scale_fill_manual(values = outcome_colors, labels = outcome_labels) +
  scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
  labs(
    title = "Figure 2 The crude rate differences of lung cancer incidence density, lung cancer mortality,\nand all-cause mortality comparing the screened with non-screened groups",
    x = NULL,
    y = "Rate Difference (per 1000 person-years)",
    fill = "Outcome"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 添加分组标签
# 手动添加分组标签在X轴上方
p <- p + annotate("text", x = 1, y = y_max * 0.95, label = "Overall", fontface = "bold", size = 3.5, hjust = 0.5) +
  annotate("text", x = 2.5, y = y_max * 0.95, label = "Sex", fontface = "bold", size = 3.5, hjust = 0.5) +
  annotate("text", x = 5, y = y_max * 0.95, label = "Age group (years)", fontface = "bold", size = 3.5, hjust = 0.5) +
  annotate("text", x = 8, y = y_max * 0.95, label = "Family history", fontface = "bold", size = 3.5, hjust = 0.5) +
  annotate("text", x = 10.5, y = y_max * 0.95, label = "Occupational hazard", fontface = "bold", size = 3.5, hjust = 0.5) +
  annotate("text", x = 13, y = y_max * 0.95, label = "Chronic respiratory disease", fontface = "bold", size = 3.5, hjust = 0.5)

# 导出为TIFF（仅TIFF，不要PNG）
tiff("/workspace/figure2_rd_barplot.tiff", width = 16, height = 10, units = "in", res = 300, compression = "lzw")
print(p)
dev.off()

cat("\nFigure 2 柱状图已导出:\n")
cat("- /workspace/figure2_rd_barplot.tiff (300 dpi, LZW压缩)\n")

cat("\n========================================\n")
cat("柱状图导出完成!\n")
cat("========================================\n")
