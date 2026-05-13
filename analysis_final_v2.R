# ============================================
# 最终版：高危筛查组 vs 高危未筛查组
# 参考图样式：左侧表格 + 右侧森林图
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
# 计算粗率、RR、95%CI和P值
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
    # P值 - 使用正态分布
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
# 定义亚组结构 - 参考图样式
# ============================================
subgroups <- list(
  list(name = "Overall", var = NULL, val = NULL, label = "Overall", level = 0),
  list(name = "Sex", var = "gender", val = "M", label = "Male", level = 1),
  list(name = "Sex", var = "gender", val = "F", label = "Female", level = 1),
  list(name = "Age group (years)", var = "age_group", val = 1, label = "40-49", level = 1),
  list(name = "Age group (years)", var = "age_group", val = 2, label = "50-54", level = 1),
  list(name = "Age group (years)", var = "age_group", val = 3, label = "55-74", level = 1),
  list(name = "Family history", var = "yijijiazushi", val = 0, label = "No", level = 1),
  list(name = "Family history", var = "yijijiazushi", val = 1, label = "Yes", level = 1),
  list(name = "Occupational hazard", var = "as_occupational_hazard_exposure", val = 0, label = "No", level = 1),
  list(name = "Occupational hazard", var = "as_occupational_hazard_exposure", val = 1, label = "Yes", level = 1),
  list(name = "Chronic respiratory disease", var = "as_chronic_respiratory_disease", val = 0, label = "No", level = 1),
  list(name = "Chronic respiratory disease", var = "as_chronic_respiratory_disease", val = 1, label = "Yes", level = 1)
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
    
    res$outcome <- outcome
    res$group <- group_name
    res$subgroup <- label
    res$level <- sg$level
    
    rr_str <- ifelse(is.na(res$rr) || is.infinite(res$rr), "N/A", sprintf("%.2f", res$rr))
    ci_str <- ifelse(is.na(res$ci_lower) || is.na(res$ci_upper), "N/A", 
                     sprintf("(%.2f-%.2f)", res$ci_lower, res$ci_upper))
    p_str <- ifelse(is.na(res$p_value), "N/A", ifelse(res$p_value < 0.0001, "<0.0001", sprintf("%.4f", res$p_value)))
    
    cat(sprintf("%s: RR=%s %s, P=%s\n", label, rr_str, ci_str, p_str))
    
    res_list[[length(res_list) + 1]] <- res
  }
  
  all_results[[outcome]] <- res_list
}

# ============================================
# 创建参考图样式的森林图
# 左侧：Group标签 + RR(95%CI) + P值
# 右侧：森林图
# ============================================

create_forest_plot <- function(dt_list, outcome_name) {
  # 准备数据框
  n <- length(dt_list)
  plot_df <- data.frame(
    y_pos = n:1,  # 反转，Overall在最上面
    group = character(n),
    subgroup = character(n),
    level = numeric(n),
    rr = numeric(n),
    ci_lower = numeric(n),
    ci_upper = numeric(n),
    p_value = numeric(n),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n) {
    plot_df$group[i] <- dt_list[[i]]$group
    plot_df$subgroup[i] <- dt_list[[i]]$subgroup
    plot_df$level[i] <- dt_list[[i]]$level
    plot_df$rr[i] <- dt_list[[i]]$rr
    plot_df$ci_lower[i] <- dt_list[[i]]$ci_lower
    plot_df$ci_upper[i] <- dt_list[[i]]$ci_upper
    plot_df$p_value[i] <- dt_list[[i]]$p_value
  }
  
  # 创建显示标签
  plot_df$display_label <- ""
  for (i in 1:n) {
    if (plot_df$subgroup[i] == "Overall") {
      plot_df$display_label[i] <- "Overall"
    } else if (plot_df$level[i] == 0) {
      plot_df$display_label[i] <- plot_df$group[i]
    } else {
      # 检查是否是该组的第一个
      if (i > 1 && plot_df$group[i] != plot_df$group[i-1] && plot_df$level[i] == 1) {
        plot_df$display_label[i] <- paste0("  ", plot_df$subgroup[i])
      } else {
        plot_df$display_label[i] <- paste0("  ", plot_df$subgroup[i])
      }
    }
  }
  
  # 创建RR文本
  plot_df$rr_text <- ""
  plot_df$ci_text <- ""
  plot_df$p_text <- ""
  
  for (i in 1:n) {
    if (is.na(plot_df$rr[i]) || is.infinite(plot_df$rr[i])) {
      plot_df$rr_text[i] <- "N/A"
      plot_df$ci_text[i] <- ""
    } else {
      plot_df$rr_text[i] <- sprintf("%.2f", plot_df$rr[i])
      if (!is.na(plot_df$ci_lower[i]) && !is.na(plot_df$ci_upper[i])) {
        plot_df$ci_text[i] <- sprintf("(%.2f-%.2f)", plot_df$ci_lower[i], plot_df$ci_upper[i])
      } else {
        plot_df$ci_text[i] <- ""
      }
    }
    
    if (is.na(plot_df$p_value[i])) {
      plot_df$p_text[i] <- "N/A"
    } else if (plot_df$p_value[i] < 0.0001) {
      plot_df$p_text[i] <- "<0.0001"
    } else {
      plot_df$p_text[i] <- sprintf("%.4f", plot_df$p_value[i])
    }
  }
  
  # 过滤有效数据用于绘图
  plot_df_valid <- plot_df[!is.na(plot_df$rr) & !is.infinite(plot_df$rr), ]
  
  color <- outcome_colors[outcome_name]
  
  # 创建森林图（右侧）
  p_forest <- ggplot() +
    # 置信区间
    geom_segment(data = plot_df_valid,
                 aes(x = ci_lower, xend = ci_upper, y = y_pos, yend = y_pos),
                 color = color, linewidth = 1.2, alpha = 0.7) +
    # RR点
    geom_point(data = plot_df_valid,
               aes(x = rr, y = y_pos),
               color = color, size = 4, shape = 21, fill = color, stroke = 0.8) +
    # 参考线 RR=1
    geom_vline(xintercept = 1, color = "black", linewidth = 0.8) +
    # 对数刻度
    scale_x_log10(limits = c(0.15, 5),
                  breaks = c(0.5, 1, 1.5, 2, 2.5),
                  labels = c("0.5", "1", "1.5", "2", "2.5")) +
    scale_y_continuous(breaks = plot_df$y_pos, labels = plot_df$display_label) +
    coord_cartesian(ylim = c(0.5, n + 0.5)) +
    labs(title = outcome_labels[outcome_name], x = "RR (95% CI)", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(color = color, size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_blank(),  # 隐藏Y轴标签，用表格替代
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(list(
    forest = p_forest,
    plot_df = plot_df,
    color = color
  ))
}

# 创建三个结局的图
res1 <- create_forest_plot(all_results[["LC_event"]], "LC_event")
res2 <- create_forest_plot(all_results[["LC_death"]], "LC_death")
res3 <- create_forest_plot(all_results[["all_cause_death"]], "all_cause_death")

# ============================================
# 创建左侧表格（使用grid/textGrob）
# ============================================

create_table_grob <- function(plot_df, color) {
  n <- nrow(plot_df)
  
  # 创建表格内容
  table_data <- data.frame(
    Group = plot_df$display_label,
    RR = plot_df$rr_text,
    CI = plot_df$ci_text,
    P = plot_df$p_text,
    stringsAsFactors = FALSE
  )
  
  # 使用tableGrob创建表格
  tg <- tableGrob(
    table_data,
    rows = NULL,
    cols = c("Group", "RR", "95% CI", "P"),
    theme = ttheme_minimal(
      core = list(
        fg_params = list(fontsize = 9, fontface = "plain", col = "black"),
        bg_params = list(fill = "white", col = NA)
      ),
      colhead = list(
        fg_params = list(fontsize = 10, fontface = "bold", col = color),
        bg_params = list(fill = "white", col = NA)
      )
    )
  )
  
  # 调整列宽
  tg$widths <- unit(c(0.4, 0.15, 0.25, 0.2), "npc")
  
  return(tg)
}

# 创建表格
tg1 <- create_table_grob(res1$plot_df, res1$color)
tg2 <- create_table_grob(res2$plot_df, res2$color)
tg3 <- create_table_grob(res3$plot_df, res3$color)

# ============================================
# 组合表格和森林图
# ============================================

# 为每个结局创建组合图
combine_plot <- function(tg, p_forest, title_text) {
  # 创建布局
  layout <- grid.layout(nrow = 1, ncol = 2, 
                        widths = unit(c(0.45, 0.55), "npc"),
                        heights = unit(1, "npc"))
  
  # 创建grob
  combined_grob <- gTree(children = gList(
    rectGrob(gp = gpar(fill = "white", col = NA)),
    tg,
    ggplotGrob(p_forest)
  ))
  
  # 使用arrangeGrob
  pg <- arrangeGrob(
    tg, ggplotGrob(p_forest),
    ncol = 2, widths = c(0.45, 0.55),
    top = textGrob(title_text, gp = gpar(fontsize = 14, fontface = "bold", col = outcome_colors[title_text]))
  )
  
  return(pg)
}

# 由于直接组合比较复杂，我们使用另一种方式：
# 在ggplot中使用annotation_custom来添加表格

create_combined_plot <- function(res, outcome_name) {
  plot_df <- res$plot_df
  color <- res$color
  
  # 过滤有效数据
  plot_df_valid <- plot_df[!is.na(plot_df$rr) & !is.infinite(plot_df$rr), ]
  
  # 创建基础图（无Y轴标签）
  p <- ggplot() +
    geom_segment(data = plot_df_valid,
                 aes(x = ci_lower, xend = ci_upper, y = y_pos, yend = y_pos),
                 color = color, linewidth = 1.5, alpha = 0.7) +
    geom_point(data = plot_df_valid,
               aes(x = rr, y = y_pos),
               color = color, size = 4, shape = 21, fill = color, stroke = 0.8) +
    geom_vline(xintercept = 1, color = "black", linewidth = 0.8) +
    scale_x_log10(limits = c(0.15, 5),
                  breaks = c(0.5, 1, 1.5, 2, 2.5),
                  labels = c("0.5", "1", "1.5", "2", "2.5")) +
    scale_y_continuous(breaks = plot_df$y_pos, 
                       labels = plot_df$display_label,
                       expand = expansion(add = c(0.5, 0.5))) +
    labs(title = outcome_labels[outcome_name], x = "RR (95% CI)", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(color = color, size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9.5, lineheight = 0.9, hjust = 0),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(10, 150, 10, 10)
    )
  
  # 添加RR数值标签
  p <- p + annotate("text", x = 3.5, y = plot_df$y_pos,
                    label = paste0(plot_df$rr_text, "\n", plot_df$ci_text),
                    hjust = 0, vjust = 0.5, size = 3, color = color, fontface = "plain")
  
  # 添加P值标签
  p <- p + annotate("text", x = 4.8, y = plot_df$y_pos,
                    label = plot_df$p_text,
                    hjust = 0, vjust = 0.5, size = 3, color = "black", fontface = "plain")
  
  return(p)
}

p1 <- create_combined_plot(res1, "LC_event")
p2 <- create_combined_plot(res2, "LC_death")
p3 <- create_combined_plot(res3, "all_cause_death")

# ============================================
# 添加列标题注释
# ============================================

# 创建带列标题的组合图
final_combined <- grid.arrange(
  p1, p2, p3,
  ncol = 3,
  top = textGrob(
    "Figure 1 The crude rate ratios of lung cancer incidence density, lung cancer mortality,\nand all-cause mortality comparing the screened with non-screened groups",
    gp = gpar(fontsize = 15, fontface = "bold")
  )
)

# 保存为TIFF
tiff("/workspace/forest_plot_final_v2.tiff", width = 28, height = 14, units = "in", res = 300, compression = "lzw")
grid.draw(final_combined)
dev.off()

png("/workspace/forest_plot_final_v2.png", width = 2800, height = 1400, res = 150)
grid.draw(final_combined)
dev.off()

cat("\n图片已保存:\n")
cat("- /workspace/forest_plot_final_v2.tiff (300 dpi, LZW压缩)\n")
cat("- /workspace/forest_plot_final_v2.png (预览)\n")

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
      CI_Upper = ifelse(is.na(res$ci_upper), NA, round(res$ci_upper, 2)),
      P_Value = ifelse(is.na(res$p_value), NA, res$p_value)
    ))
  }
}

fwrite(result_table, "/workspace/rr_results_final_v2.csv")
cat("结果表格已保存到 /workspace/rr_results_final_v2.csv\n")
