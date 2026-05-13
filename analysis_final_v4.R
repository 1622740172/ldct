# ============================================
# 最终版：高危筛查组 vs 高危未筛查组
# 参考图样式：左侧分组 + 右侧森林图 + RR/CI/P值
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
# 计算RR
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
    z <- log_rr / se_log_rr
    p_value <- 2 * (1 - pnorm(abs(z)))
  } else {
    ci_lower <- NA
    ci_upper <- NA
    p_value <- NA
  }
  
  return(list(
    screened_events = e1, unscreened_events = e0,
    screened_py = py1, unscreened_py = py0,
    screened_rate = rate1, unscreened_rate = rate0,
    rr = rr, ci_lower = ci_lower, ci_upper = ci_upper,
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
        rr = NA, ci_lower = NA, ci_upper = NA, p_value = NA
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
      
      res <- calculate_rr(dt_sub, outcome, time_col)
      res$type <- "data"
      res$group <- group_name
      res$subgroup <- label
      
      rr_str <- ifelse(is.na(res$rr) || is.infinite(res$rr), "N/A", sprintf("%.2f", res$rr))
      ci_str <- ifelse(is.na(res$ci_lower) || is.na(res$ci_upper), "N/A", 
                       sprintf("(%.2f-%.2f)", res$ci_lower, res$ci_upper))
      p_str <- ifelse(is.na(res$p_value), "N/A", ifelse(res$p_value < 0.0001, "<0.0001", sprintf("%.4f", res$p_value)))
      
      cat(sprintf("%s: RR=%s %s, P=%s\n", label, rr_str, ci_str, p_str))
    }
    
    res_list[[length(res_list) + 1]] <- res
  }
  
  all_results[[outcome]] <- res_list
}

# ============================================
# 创建森林图
# ============================================

create_forest_plot <- function(dt_list, outcome_name) {
  n <- length(dt_list)
  
  plot_df <- data.frame(
    y_pos = n:1,
    type = character(n),
    group = character(n),
    subgroup = character(n),
    rr = numeric(n),
    ci_lower = numeric(n),
    ci_upper = numeric(n),
    p_value = numeric(n),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n) {
    plot_df$type[i] <- dt_list[[i]]$type
    plot_df$group[i] <- dt_list[[i]]$group
    plot_df$subgroup[i] <- dt_list[[i]]$subgroup
    plot_df$rr[i] <- ifelse(is.null(dt_list[[i]]$rr), NA, dt_list[[i]]$rr)
    plot_df$ci_lower[i] <- ifelse(is.null(dt_list[[i]]$ci_lower), NA, dt_list[[i]]$ci_lower)
    plot_df$ci_upper[i] <- ifelse(is.null(dt_list[[i]]$ci_upper), NA, dt_list[[i]]$ci_upper)
    plot_df$p_value[i] <- ifelse(is.null(dt_list[[i]]$p_value), NA, dt_list[[i]]$p_value)
  }
  
  # 显示标签
  plot_df$display_label <- ""
  for (i in 1:n) {
    if (plot_df$type[i] == "header") {
      plot_df$display_label[i] <- plot_df$group[i]
    } else if (plot_df$subgroup[i] == "Overall") {
      plot_df$display_label[i] <- "Overall"
    } else {
      plot_df$display_label[i] <- paste0("  ", plot_df$subgroup[i])
    }
  }
  
  # 文本
  plot_df$rr_ci_text <- ""
  plot_df$p_text <- ""
  
  for (i in 1:n) {
    if (plot_df$type[i] == "header") {
      plot_df$rr_ci_text[i] <- ""
      plot_df$p_text[i] <- ""
    } else if (is.na(plot_df$rr[i]) || is.infinite(plot_df$rr[i])) {
      plot_df$rr_ci_text[i] <- "N/A"
      plot_df$p_text[i] <- "N/A"
    } else {
      rr_str <- sprintf("%.2f", plot_df$rr[i])
      ci_str <- sprintf("(%.2f-%.2f)", plot_df$ci_lower[i], plot_df$ci_upper[i])
      plot_df$rr_ci_text[i] <- paste0(rr_str, " ", ci_str)
      
      if (plot_df$p_value[i] < 0.0001) {
        plot_df$p_text[i] <- "<0.0001"
      } else {
        plot_df$p_text[i] <- sprintf("%.4f", plot_df$p_value[i])
      }
    }
  }
  
  color <- outcome_colors[outcome_name]
  plot_df_valid <- plot_df[plot_df$type == "data" & !is.na(plot_df$rr) & !is.infinite(plot_df$rr), ]
  
  # 计算需要显示的x范围
  all_ci <- c(plot_df_valid$ci_lower, plot_df_valid$ci_upper)
  all_ci <- all_ci[!is.na(all_ci) & !is.infinite(all_ci)]
  x_min <- min(0.2, min(all_ci) * 0.8)
  x_max <- max(5, max(all_ci) * 1.2)
  
  # 创建森林图
  p <- ggplot() +
    geom_segment(data = plot_df_valid,
                 aes(x = ci_lower, xend = ci_upper, y = y_pos, yend = y_pos),
                 color = color, linewidth = 1.5, alpha = 0.7) +
    geom_point(data = plot_df_valid,
               aes(x = rr, y = y_pos),
               color = color, size = 4, shape = 21, fill = color, stroke = 0.8) +
    geom_vline(xintercept = 1, color = "black", linewidth = 0.8) +
    scale_x_log10(limits = c(x_min, x_max),
                  breaks = c(0.5, 1, 1.5, 2, 2.5),
                  labels = c("0.5", "1", "1.5", "2", "2.5")) +
    scale_y_continuous(breaks = plot_df$y_pos, 
                       labels = plot_df$display_label,
                       expand = expansion(add = c(0.5, 0.5))) +
    labs(title = outcome_labels[outcome_name], x = "RR (95% CI)", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(color = color, size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10, lineheight = 0.9, hjust = 0,
                                  face = ifelse(plot_df$type == "header", "bold", "plain")),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(10, 250, 10, 150)
    )
  
  # 标注RR+CI
  p <- p + annotate("text", x = x_max * 1.15, y = plot_df$y_pos,
                    label = plot_df$rr_ci_text,
                    hjust = 0, vjust = 0.5, size = 3.2, 
                    color = ifelse(plot_df$type == "header", "white", color),
                    fontface = ifelse(plot_df$type == "header", "bold", "plain"))
  
  # 标注P值
  p <- p + annotate("text", x = x_max * 1.65, y = plot_df$y_pos,
                    label = plot_df$p_text,
                    hjust = 0, vjust = 0.5, size = 3.2, 
                    color = ifelse(plot_df$type == "header", "white", "black"),
                    fontface = ifelse(plot_df$type == "header", "bold", "plain"))
  
  # 添加列标题
  p <- p + annotate("text", x = x_max * 1.15, y = max(plot_df$y_pos) + 1.5,
                    label = "RR (95% CI)",
                    hjust = 0, vjust = 0.5, size = 3.5, fontface = "bold", color = color) +
    annotate("text", x = x_max * 1.65, y = max(plot_df$y_pos) + 1.5,
             label = "P value",
             hjust = 0, vjust = 0.5, size = 3.5, fontface = "bold", color = color)
  
  return(p)
}

# 创建三个结局的图
p1 <- create_forest_plot(all_results[["LC_event"]], "LC_event")
p2 <- create_forest_plot(all_results[["LC_death"]], "LC_death")
p3 <- create_forest_plot(all_results[["all_cause_death"]], "all_cause_death")

# 组合图
final_combined <- grid.arrange(
  p1, p2, p3,
  ncol = 3,
  top = textGrob(
    "Figure 1 The crude rate ratios of lung cancer incidence density, lung cancer mortality,\nand all-cause mortality comparing the screened with non-screened groups",
    gp = gpar(fontsize = 15, fontface = "bold")
  )
)

# 保存为TIFF
tiff("/workspace/forest_plot_final_v4.tiff", width = 30, height = 16, units = "in", res = 300, compression = "lzw")
grid.draw(final_combined)
dev.off()

png("/workspace/forest_plot_final_v4.png", width = 3000, height = 1600, res = 150)
grid.draw(final_combined)
dev.off()

cat("\n图片已保存:\n")
cat("- /workspace/forest_plot_final_v4.tiff (300 dpi, LZW压缩)\n")
cat("- /workspace/forest_plot_final_v4.png (预览)\n")

# 导出结果表格
result_table <- data.frame()
for (outcome in names(outcomes)) {
  for (res in all_results[[outcome]]) {
    if (res$type == "data") {
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
        CI_Upper = ifelse(is.na(res$ci_upper), NA, round(res$ci_upper, 2)),
        P_Value = ifelse(is.na(res$p_value), NA, res$p_value)
      ))
    }
  }
}

fwrite(result_table, "/workspace/rr_results_final_v4.csv")
cat("结果表格已保存到 /workspace/rr_results_final_v4.csv\n")
