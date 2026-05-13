# ============================================
# 方案B：以低危组为参照，三组比较
# ============================================

library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)

# 读取数据
df <- fread("/workspace/匹配结果_最终修正版.csv", encoding="UTF-8")

# 计算人年
df[, cancer_time_years := cancer_time / 365.25]
df[, death_time_years := death_time / 365.25]

cat("总数据:", nrow(df), "\n")

# ============================================
# 计算粗率和RR
# ============================================
calculate_rr <- function(dt, outcome_col, time_col, target_val, ref_val) {
  target <- dt[as_lc_risk == target_val]
  ref <- dt[as_lc_risk == ref_val]
  
  e_target <- target[, sum(get(outcome_col))]
  e_ref <- ref[, sum(get(outcome_col))]
  py_target <- target[, sum(get(time_col))]
  py_ref <- ref[, sum(get(time_col))]
  
  rate_target <- ifelse(py_target > 0, (e_target / py_target) * 1000, 0)
  rate_ref <- ifelse(py_ref > 0, (e_ref / py_ref) * 1000, 0)
  rr <- ifelse(rate_ref > 0, rate_target / rate_ref, NA)
  
  if (e_target > 0 && e_ref > 0 && !is.na(rr) && !is.infinite(rr)) {
    log_rr <- log(rr)
    se_log_rr <- sqrt(1/e_target + 1/e_ref)
    ci_lower <- exp(log_rr - 1.96 * se_log_rr)
    ci_upper <- exp(log_rr + 1.96 * se_log_rr)
  } else {
    ci_lower <- NA
    ci_upper <- NA
  }
  
  return(list(
    target_events = e_target, ref_events = e_ref,
    target_py = py_target, ref_py = py_ref,
    target_rate = rate_target, ref_rate = rate_ref,
    rr = rr, ci_lower = ci_lower, ci_upper = ci_upper
  ))
}

# ============================================
# 定义三组比较
# ============================================
comparisons <- list(
  list(name = "High-risk screened vs Low-risk", target = 1, ref = 0, color = "#2E86AB"),
  list(name = "High-risk non-screened vs Low-risk", target = 2, ref = 0, color = "#A23B72"),
  list(name = "High-risk screened vs High-risk non-screened", target = 1, ref = 2, color = "#F18F01")
)

# ============================================
# 定义亚组
# ============================================
subgroups <- list(
  list(name = "Overall", var = NULL, val = NULL, label = "Overall", is_header = FALSE),
  list(name = "Sex", var = "gender", val = "M", label = "Male", is_header = TRUE),
  list(name = "Sex", var = "gender", val = "F", label = "Female", is_header = FALSE),
  list(name = "Age group", var = "age_group", val = 1, label = "40-49", is_header = TRUE),
  list(name = "Age group", var = "age_group", val = 2, label = "50-54", is_header = FALSE),
  list(name = "Age group", var = "age_group", val = 3, label = "55-74", is_header = FALSE),
  list(name = "Family history", var = "yijijiazushi", val = 0, label = "No", is_header = TRUE),
  list(name = "Family history", var = "yijijiazushi", val = 1, label = "Yes", is_header = FALSE),
  list(name = "Occupational hazard", var = "as_occupational_hazard_exposure", val = 0, label = "No", is_header = TRUE),
  list(name = "Occupational hazard", var = "as_occupational_hazard_exposure", val = 1, label = "Yes", is_header = FALSE),
  list(name = "Chronic resp. disease", var = "as_chronic_respiratory_disease", val = 0, label = "No", is_header = TRUE),
  list(name = "Chronic resp. disease", var = "as_chronic_respiratory_disease", val = 1, label = "Yes", is_header = FALSE)
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

comp_colors <- c("#2E86AB", "#A23B72", "#F18F01")
comp_names <- c("High-risk screened vs Low-risk", 
                "High-risk non-screened vs Low-risk",
                "High-risk screened vs High-risk non-screened")

# ============================================
# 计算所有结果
# ============================================
all_results <- list()

for (outcome in names(outcomes)) {
  time_col <- outcomes[[outcome]]
  res_list <- list()
  
  for (sg in subgroups) {
    if (is.null(sg$var)) {
      dt_sub <- df
      label <- "Overall"
    } else {
      dt_sub <- df[get(sg$var) == sg$val]
      label <- sg$label
    }
    
    for (comp in comparisons) {
      res <- calculate_rr(dt_sub, outcome, time_col, comp$target, comp$ref)
      res$outcome <- outcome
      res$group <- sg$name
      res$subgroup <- label
      res$is_header <- sg$is_header
      res$comparison <- comp$name
      res$comp_idx <- which(sapply(comparisons, function(x) x$name == comp$name))
      res_list[[length(res_list) + 1]] <- res
    }
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
      comparison = res$comparison,
      comp_idx = res$comp_idx,
      rr = res$rr,
      ci_lower = res$ci_lower,
      ci_upper = res$ci_upper
    ))
  }
}

plot_data <- as.data.table(plot_data)

# 排序
plot_data[, group_order := 0]
plot_data[group == "Overall", group_order := 1]
plot_data[group == "Sex", group_order := 2]
plot_data[group == "Age group", group_order := 3]
plot_data[group == "Family history", group_order := 4]
plot_data[group == "Occupational hazard", group_order := 5]
plot_data[group == "Chronic resp. disease", group_order := 6]

plot_data[, subgroup_order := 0]
plot_data[subgroup == "Overall", subgroup_order := 1]
plot_data[subgroup == "Male", subgroup_order := 1]
plot_data[subgroup == "Female", subgroup_order := 2]
plot_data[subgroup == "40-49", subgroup_order := 1]
plot_data[subgroup == "50-54", subgroup_order := 2]
plot_data[subgroup == "55-74", subgroup_order := 3]
plot_data[subgroup == "No", subgroup_order := 1]
plot_data[subgroup == "Yes", subgroup_order := 2]

setorder(plot_data, outcome, group_order, subgroup_order, comp_idx)

plot_data[, grp_id := .GRP, by = .(outcome, group, subgroup)]
plot_data[, y_pos := grp_id, by = outcome]
max_y <- plot_data[, max(y_pos)]
plot_data[, y_pos := max_y - y_pos + 1]

plot_data[, display_label := ""]
plot_data[subgroup == "Overall", display_label := "Overall"]
plot_data[is_header == TRUE & subgroup != "Overall", display_label := group]
plot_data[is_header == FALSE, display_label := paste0("    ", subgroup)]

plot_data[, x_offset := (comp_idx - 2) * 0.12]

# ============================================
# 创建森林图
# ============================================
create_forest_plot <- function(dt, outcome_name) {
  dt_sub <- dt[outcome == outcome_name]
  
  dt_sub[, rr_text := ifelse(is.na(rr) | is.infinite(rr), "N/A", sprintf("%.2f", rr))]
  dt_sub[, ci_text := ifelse(is.na(ci_lower) | is.na(ci_upper), "N/A", 
                              sprintf("(%.2f-%.2f)", ci_lower, ci_upper))]
  dt_sub[, full_text := paste0(rr_text, "\n", ci_text)]
  
  dt_plot <- dt_sub[!is.na(rr) & !is.infinite(rr)]
  
  label_df <- dt_sub[, .(display_label = display_label[1]), by = y_pos]
  label_df <- label_df[order(y_pos)]
  
  p <- ggplot() +
    geom_segment(data = dt_plot, 
                 aes(x = ci_lower, xend = ci_upper, y = y_pos + x_offset, yend = y_pos + x_offset, color = comparison),
                 linewidth = 1.2, alpha = 0.7) +
    geom_point(data = dt_plot,
               aes(x = rr, y = y_pos + x_offset, color = comparison),
               size = 3, shape = 21, fill = "white", stroke = 1.2) +
    geom_vline(xintercept = 1, color = "black", linewidth = 0.8) +
    scale_x_log10(limits = c(0.1, 10),
                  breaks = c(0.5, 1, 2, 5),
                  labels = c("0.5", "1", "2", "5")) +
    scale_y_continuous(breaks = label_df$y_pos, labels = label_df$display_label) +
    coord_cartesian(ylim = c(min(label_df$y_pos) - 0.5, max(label_df$y_pos) + 0.5)) +
    labs(title = outcome_labels[outcome_name], x = "RR (95% CI)", y = NULL) +
    scale_color_manual(values = setNames(comp_colors, comp_names),
                       name = "Comparison",
                       labels = c("HR screened vs Low-risk", 
                                  "HR non-screened vs Low-risk",
                                  "HR screened vs HR non-screened")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9, lineheight = 0.85),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(10, 120, 10, 10)
    )
  
  p <- p + annotate("text", x = 6, y = dt_sub$y_pos + dt_sub$x_offset, 
                    label = dt_sub$full_text, 
                    hjust = 0, vjust = 0.5, size = 2.8, 
                    color = comp_colors[dt_sub$comp_idx], fontface = "plain")
  
  return(p)
}

p1 <- create_forest_plot(plot_data, "LC_event")
p2 <- create_forest_plot(plot_data, "LC_death")
p3 <- create_forest_plot(plot_data, "all_cause_death")

# 组合图
combined <- grid.arrange(
  p1, p2, p3,
  ncol = 3,
  top = textGrob(
    "Figure 1 The crude rate ratios of lung cancer incidence density, lung cancer mortality,\nand all-cause mortality comparing the screened with non-screened groups",
    gp = gpar(fontsize = 15, fontface = "bold")
  )
)

# 保存为TIFF
tiff("/workspace/forest_plot_scheme_B.tiff", width = 28, height = 18, units = "in", res = 300, compression = "lzw")
grid.draw(combined)
dev.off()

png("/workspace/forest_plot_scheme_B.png", width = 2800, height = 1800, res = 150)
grid.draw(combined)
dev.off()

cat("\n图片已保存:\n")
cat("- /workspace/forest_plot_scheme_B.tiff (300 dpi, LZW压缩)\n")
cat("- /workspace/forest_plot_scheme_B.png (预览)\n")

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
      Comparison = res$comparison,
      Target_Events = res$target_events,
      Target_PY = round(res$target_py, 1),
      Target_Rate = round(res$target_rate, 2),
      Ref_Events = res$ref_events,
      Ref_PY = round(res$ref_py, 1),
      Ref_Rate = round(res$ref_rate, 2),
      RR = ifelse(is.na(res$rr) || is.infinite(res$rr), NA, round(res$rr, 2)),
      CI_Lower = ifelse(is.na(res$ci_lower), NA, round(res$ci_lower, 2)),
      CI_Upper = ifelse(is.na(res$ci_upper), NA, round(res$ci_upper, 2))
    ))
  }
}

fwrite(result_table, "/workspace/rr_results_scheme_B.csv")
cat("结果表格已保存到 /workspace/rr_results_scheme_B.csv\n")
