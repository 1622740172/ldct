# ============================================
# Figure 2: Rate Difference (RD) 森林图
# The crude rate differences of lung cancer incidence density,
# lung cancer mortality, and all-cause mortality
# comparing the screened with non-screened groups
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

cat("高危人群总数:", nrow(df_high), "\n")

# ============================================
# 计算Rate Difference和详细统计量
# RD = Rate_screened - Rate_unscreened (per 1000 person-years)
# SE_RD = sqrt(e1/py1^2 + e0/py0^2) * 1000
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
    z <- rd / se_rd
    p_value <- 2 * (1 - pnorm(abs(z)))
  } else {
    se_rd <- NA
    ci_lower <- NA
    ci_upper <- NA
    p_value <- NA
  }
  
  return(list(
    n_screened = n1, n_unscreened = n0,
    events_screened = e1, events_unscreened = e0,
    py_screened = py1, py_unscreened = py0,
    rate_screened = rate1, rate_unscreened = rate0,
    rd = rd, ci_lower = ci_lower, ci_upper = ci_upper,
    p_value = p_value
  ))
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
  LC_event = "Lung cancer incidence",
  LC_death = "Lung cancer mortality", 
  all_cause_death = "All-cause mortality"
)

outcome_colors <- c(LC_event = "#2E86AB", LC_death = "#A23B72", all_cause_death = "#F18F01")

# ============================================
# 计算所有结果
# ============================================
all_results <- list()

for (outcome in names(outcomes)) {
  time_col <- outcomes[[outcome]]
  cat("\n========================================\n")
  cat("Outcome:", outcome, "\n")
  cat("========================================\n")
  
  res_list <- list()
  
  for (sg in subgroups) {
    if (sg$type == "header") {
      res <- list(
        type = "header", group = sg$name, subgroup = sg$name,
        n_screened = NA, n_unscreened = NA,
        events_screened = NA, events_unscreened = NA,
        rd = NA, ci_lower = NA, ci_upper = NA, p_value = NA
      )
    } else {
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
      res$type <- "data"
      res$group <- group_name
      res$subgroup <- label
      
      rd_str <- ifelse(is.na(res$rd), "N/A", sprintf("%.2f", res$rd))
      ci_str <- ifelse(is.na(res$ci_lower) || is.na(res$ci_upper), "N/A", 
                       sprintf("(%.2f-%.2f)", res$ci_lower, res$ci_upper))
      p_str <- ifelse(is.na(res$p_value), "N/A", ifelse(res$p_value < 0.0001, "<0.0001", sprintf("%.4f", res$p_value)))
      
      cat(sprintf("%s: n_screened=%d, n_unscreened=%d, RD=%s %s, P=%s\n", 
                  label, res$n_screened, res$n_unscreened, rd_str, ci_str, p_str))
    }
    
    res_list[[length(res_list) + 1]] <- res
  }
  
  all_results[[outcome]] <- res_list
}

# ============================================
# 创建数据框
# ============================================
n <- length(subgroups)
plot_df <- data.frame(
  y_pos = n:1,
  type = character(n),
  group = character(n),
  subgroup = character(n),
  display_label = character(n),
  stringsAsFactors = FALSE
)

for (i in 1:n) {
  plot_df$type[i] <- subgroups[[i]]$type
  plot_df$group[i] <- subgroups[[i]]$name
  if (subgroups[[i]]$type == "header") {
    plot_df$subgroup[i] <- subgroups[[i]]$name
    plot_df$display_label[i] <- subgroups[[i]]$name
  } else {
    plot_df$subgroup[i] <- subgroups[[i]]$label
    if (subgroups[[i]]$label == "Overall") {
      plot_df$display_label[i] <- "Overall"
    } else {
      plot_df$display_label[i] <- paste0("  ", subgroups[[i]]$label)
    }
  }
}

# 为每个结局准备数据
for (outcome in names(outcomes)) {
  res_list <- all_results[[outcome]]
  
  rd_col <- paste0("rd_", outcome)
  ci_lower_col <- paste0("ci_lower_", outcome)
  ci_upper_col <- paste0("ci_upper_", outcome)
  p_col <- paste0("p_", outcome)
  
  plot_df[[rd_col]] <- NA
  plot_df[[ci_lower_col]] <- NA
  plot_df[[ci_upper_col]] <- NA
  plot_df[[p_col]] <- NA
  
  for (i in 1:n) {
    res <- res_list[[i]]
    if (res$type == "data") {
      plot_df[[rd_col]][i] <- res$rd
      plot_df[[ci_lower_col]][i] <- res$ci_lower
      plot_df[[ci_upper_col]][i] <- res$ci_upper
      plot_df[[p_col]][i] <- res$p_value
    }
  }
}

# ============================================
# 创建RD森林图 - 左侧只显示一次变量名
# ============================================

create_rd_forest_plot <- function(plot_df, outcome_name, show_y_labels = TRUE, show_rd_text = TRUE) {
  color <- outcome_colors[outcome_name]
  
  rd_col <- paste0("rd_", outcome_name)
  ci_lower_col <- paste0("ci_lower_", outcome_name)
  ci_upper_col <- paste0("ci_upper_", outcome_name)
  
  valid_idx <- !is.na(plot_df[[rd_col]])
  plot_df_valid <- plot_df[valid_idx, ]
  
  all_ci <- c(plot_df[[ci_lower_col]], plot_df[[ci_upper_col]])
  all_ci <- all_ci[!is.na(all_ci)]
  x_min <- min(-5, min(all_ci) * 1.2)
  x_max <- max(5, max(all_ci) * 1.2)
  
  # 准备RD文本标签
  plot_df$rd_text <- ""
  for (i in 1:nrow(plot_df)) {
    if (plot_df$type[i] == "data" && !is.na(plot_df[[rd_col]][i])) {
      rd_val <- plot_df[[rd_col]][i]
      ci_l <- plot_df[[ci_lower_col]][i]
      ci_u <- plot_df[[ci_upper_col]][i]
      plot_df$rd_text[i] <- sprintf("%.2f (%.2f-%.2f)", rd_val, ci_l, ci_u)
    }
  }
  
  p <- ggplot() +
    geom_segment(data = plot_df_valid,
                 aes(x = .data[[ci_lower_col]], xend = .data[[ci_upper_col]], y = y_pos, yend = y_pos),
                 color = color, linewidth = 1.5, alpha = 0.7) +
    geom_point(data = plot_df_valid,
               aes(x = .data[[rd_col]], y = y_pos),
               color = color, size = 4, shape = 21, fill = color, stroke = 0.8) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
    scale_x_continuous(limits = c(x_min, x_max)) +
    labs(title = outcome_labels[outcome_name], x = "RD (95% CI)", y = NULL) +
    theme_minimal()
  
  if (show_y_labels) {
    p <- p + scale_y_continuous(breaks = plot_df$y_pos, 
                                labels = plot_df$display_label,
                                expand = expansion(add = c(0.5, 0.5)))
  } else {
    p <- p + scale_y_continuous(breaks = plot_df$y_pos, 
                                labels = NULL,
                                expand = expansion(add = c(0.5, 0.5)))
  }
  
  # 添加RD数值标签在右侧
  if (show_rd_text) {
    p <- p + geom_text(data = plot_df[plot_df$type == "data", ],
                       aes(x = x_max * 0.98, y = y_pos, label = rd_text),
                       hjust = 1, vjust = 0.5, size = 2.8, color = "black")
  }
  
  p <- p + theme(
    plot.title = element_text(color = color, size = 13, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10, lineheight = 0.9, hjust = 0,
                                face = ifelse(plot_df$type == "header", "bold", "plain")),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(10, 60, 10, 10)
  )
  
  return(p)
}

# 创建三个RD森林图
p1 <- create_rd_forest_plot(plot_df, "LC_event", show_y_labels = TRUE, show_rd_text = TRUE)
p2 <- create_rd_forest_plot(plot_df, "LC_death", show_y_labels = FALSE, show_rd_text = TRUE)
p3 <- create_rd_forest_plot(plot_df, "all_cause_death", show_y_labels = FALSE, show_rd_text = TRUE)

# ============================================
# 使用cowplot组合（避免grid.arrange的遮盖问题）
# ============================================

tryCatch({
  library(cowplot)
  
  # 创建标题
  title_grob <- ggdraw() + 
    draw_label("Figure 2 The crude rate differences of lung cancer incidence density, lung cancer mortality,\nand all-cause mortality comparing the screened with non-screened groups",
               fontface = "bold", size = 14, hjust = 0.5)
  
  # 组合三个图，调整宽度比例给左侧标签留空间
  plot_row <- plot_grid(p1, p2, p3, 
                        ncol = 3, 
                        rel_widths = c(1.6, 1, 1),
                        align = "h", axis = "tb")
  
  forest_combined <- plot_grid(title_grob, plot_row, ncol = 1, rel_heights = c(0.12, 1))
  
  # 导出为TIFF
  tiff("/workspace/figure2_rd_forest.tiff", width = 26, height = 12, units = "in", res = 300, compression = "lzw")
  print(forest_combined)
  dev.off()
  
  png("/workspace/figure2_rd_forest.png", width = 2600, height = 1200, res = 150)
  print(forest_combined)
  dev.off()
  
  cat("\nFigure 2 森林图已导出 (使用cowplot):\n")
  cat("- /workspace/figure2_rd_forest.tiff (300 dpi, LZW压缩)\n")
  cat("- /workspace/figure2_rd_forest.png (预览)\n")
  
}, error = function(e) {
  cat("\ncowplot不可用，使用grid.arrange...\n")
  
  # 使用grid.arrange，但增加整体宽度和调整比例
  forest_combined <- grid.arrange(
    p1, p2, p3,
    ncol = 3,
    widths = c(1.6, 1, 1),
    top = textGrob(
      "Figure 2 The crude rate differences of lung cancer incidence density, lung cancer mortality,\nand all-cause mortality comparing the screened with non-screened groups",
      gp = gpar(fontsize = 14, fontface = "bold")
    )
  )
  
  tiff("/workspace/figure2_rd_forest.tiff", width = 26, height = 12, units = "in", res = 300, compression = "lzw")
  grid.draw(forest_combined)
  dev.off()
  
  png("/workspace/figure2_rd_forest.png", width = 2600, height = 1200, res = 150)
  grid.draw(forest_combined)
  dev.off()
  
  cat("\nFigure 2 森林图已导出 (使用grid.arrange):\n")
  cat("- /workspace/figure2_rd_forest.tiff (300 dpi, LZW压缩)\n")
  cat("- /workspace/figure2_rd_forest.png (预览)\n")
})

# ============================================
# 导出结果表格（CSV）
# ============================================
result_table <- data.frame()
for (outcome in names(outcomes)) {
  for (res in all_results[[outcome]]) {
    if (res$type == "data") {
      result_table <- rbind(result_table, data.frame(
        Outcome = outcome_labels[outcome],
        Group = res$group,
        Subgroup = res$subgroup,
        Screened_N = res$n_screened,
        Screened_Events = res$events_screened,
        Unscreened_N = res$n_unscreened,
        Unscreened_Events = res$events_unscreened,
        Rate_Screened = ifelse(is.na(res$rate_screened), NA, round(res$rate_screened, 2)),
        Rate_Unscreened = ifelse(is.na(res$rate_unscreened), NA, round(res$rate_unscreened, 2)),
        RD = ifelse(is.na(res$rd), NA, round(res$rd, 2)),
        CI_Lower = ifelse(is.na(res$ci_lower), NA, round(res$ci_lower, 2)),
        CI_Upper = ifelse(is.na(res$ci_upper), NA, round(res$ci_upper, 2)),
        P_Value = ifelse(is.na(res$p_value), NA, res$p_value)
      ))
    }
  }
}

fwrite(result_table, "/workspace/rd_results_figure2.csv")
cat("\n结果表格已保存到 /workspace/rd_results_figure2.csv\n")

cat("\n========================================\n")
cat("Figure 2 所有导出完成!\n")
cat("========================================\n")
