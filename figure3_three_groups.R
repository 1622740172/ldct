# ============================================
# Figure 3: Three-group cumulative curves with competing risk adjustment
# Groups: Low-risk (as_lc_risk=0), Screened (as_lc_risk=1), Non-screened (as_lc_risk=2)
# (A) Lung cancer incidence
# (B) Lung cancer mortality  
# (C) All-cause mortality
# X-axis: years, 1-year intervals
# ============================================

library(data.table)
library(ggplot2)
library(cowplot)

# 读取数据
df <- fread("/workspace/匹配结果_最终修正版.csv", encoding="UTF-8")

# 定义三组
df[, group := factor(as_lc_risk, 
                     levels = c(0, 1, 2),
                     labels = c("Low-risk", "Screened", "Non-screened"))]

# 转换时间为年
df[, time_years := cancer_time / 365.25]

# 确保事件变量是数值型
df[, LC_event := as.numeric(LC_event)]
df[, LC_death := as.numeric(LC_death)]
df[, all_cause_death := as.numeric(all_cause_death)]

# ============================================
# 竞争风险模型：Aalen-Johansen 估计器
# ============================================

calculate_aj_cumulative <- function(dt, event_col, competing_col, time_col, max_years = 12) {
  groups <- unique(dt$group)
  result <- data.frame()
  summary_table <- data.frame()
  
  for (g in groups) {
    sub_dt <- dt[group == g]
    n_total <- nrow(sub_dt)
    
    all_times <- sort(unique(sub_dt[[time_col]]))
    all_times <- all_times[all_times > 0]
    
    cum_inc <- 0
    surv_overall <- 1
    
    time_points <- c(0)
    cum_values <- c(0)
    
    for (t in all_times) {
      at_risk <- sum(sub_dt[[time_col]] >= t)
      
      if (at_risk == 0) next
      
      d_event <- sum(sub_dt[[time_col]] == t & sub_dt[[event_col]] == 1)
      d_competing <- sum(sub_dt[[time_col]] == t & sub_dt[[competing_col]] == 1)
      d_total <- d_event + d_competing
      
      if (at_risk > 0) {
        surv_overall <- surv_overall * (1 - d_total / at_risk)
        cum_inc <- cum_inc + surv_overall * (d_event / at_risk)
      }
      
      time_points <- c(time_points, t / 365.25)
      cum_values <- c(cum_values, cum_inc * 100)
    }
    
    # 插值到标准时间网格
    time_grid <- seq(0, max_years, by = 0.05)
    cum_interp <- approx(time_points, cum_values, xout = time_grid, 
                         method = "constant", f = 0, yleft = 0)$y
    cum_interp[is.na(cum_interp)] <- max(cum_interp, na.rm = TRUE)
    
    result <- rbind(result, data.frame(
      time = time_grid,
      cumulative = cum_interp,
      group = g,
      stringsAsFactors = FALSE
    ))
    
    # 计算1年、5年、10年的累积率
    for (yr in c(1, 5, 10)) {
      idx <- which.min(abs(time_grid - yr))
      summary_table <- rbind(summary_table, data.frame(
        group = g,
        year = yr,
        cumulative_rate = round(cum_interp[idx], 2),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(list(curve = result, summary = summary_table))
}

# ============================================
# 计算三个结局的累积曲线（三组+竞争风险调整）
# ============================================

cat("计算肺癌发病率累积曲线（三组+竞争风险调整）...\n")
aj_lc_event <- calculate_aj_cumulative(
  df, 
  event_col = "LC_event", 
  competing_col = "all_cause_death",
  time_col = "cancer_time",
  max_years = 12
)

cat("计算肺癌死亡率累积曲线（三组+竞争风险调整）...\n")
df[, non_lc_death := ifelse(all_cause_death == 1 & LC_death == 0, 1, 0)]
aj_lc_death <- calculate_aj_cumulative(
  df, 
  event_col = "LC_death", 
  competing_col = "non_lc_death",
  time_col = "death_time",
  max_years = 12
)

cat("计算全因死亡率累积曲线（三组）...\n")
aj_all_death <- calculate_aj_cumulative(
  df, 
  event_col = "all_cause_death", 
  competing_col = "LC_event",
  time_col = "death_time",
  max_years = 12
)

# ============================================
# 设置颜色和样式 - 三组清晰区分
# ============================================
group_colors <- c("Low-risk" = "#2E86AB", "Screened" = "#A23B72", "Non-screened" = "#F18F01")
group_linetypes <- c("Low-risk" = "solid", "Screened" = "dashed", "Non-screened" = "dotted")

# ============================================
# 创建累积曲线图函数（带清晰图例和数值标注）
# ============================================

create_cumulative_plot <- function(aj_result, title_text, y_label, y_max = NULL) {
  cum_data <- aj_result$curve
  summary_data <- aj_result$summary
  
  # 确保group是因子且顺序正确
  cum_data$group <- factor(cum_data$group, levels = c("Low-risk", "Screened", "Non-screened"))
  summary_data$group <- factor(summary_data$group, levels = c("Low-risk", "Screened", "Non-screened"))
  
  if (is.null(y_max)) {
    y_max <- max(cum_data$cumulative, na.rm = TRUE) * 1.6
  }
  
  # 准备标注文本（按组显示）
  label_lines <- c()
  for (g in levels(summary_data$group)) {
    sub <- summary_data[summary_data$group == g, ]
    rates <- paste0(sub$year, "y=", sub$cumulative_rate, "%")
    label_lines <- c(label_lines, paste0(g, ": ", paste(rates, collapse = ", ")))
  }
  label_text <- paste(label_lines, collapse = "\n")
  
  p <- ggplot(cum_data, aes(x = time, y = cumulative, color = group, linetype = group)) +
    geom_line(linewidth = 1.5) +
    scale_color_manual(
      values = group_colors,
      labels = c("Low-risk (as_lc_risk=0)", "Screened (as_lc_risk=1)", "Non-screened (as_lc_risk=2)")
    ) +
    scale_linetype_manual(
      values = group_linetypes,
      labels = c("Low-risk (as_lc_risk=0)", "Screened (as_lc_risk=1)", "Non-screened (as_lc_risk=2)")
    ) +
    scale_x_continuous(
      breaks = seq(0, 12, by = 1),
      labels = seq(0, 12, by = 1),
      limits = c(0, 12),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = c(0, 0),
      labels = function(x) paste0(x, "%")
    ) +
    labs(
      title = title_text,
      x = "Follow-up time (years)",
      y = y_label,
      color = "Group",
      linetype = "Group"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
      legend.box.background = element_rect(color = "gray80", linewidth = 0.5),
      legend.key.size = unit(1.2, "cm"),
      legend.key.width = unit(1.5, "cm"),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  # 添加数值标注在图右下角
  p <- p + annotate("text", x = 7.5, y = y_max * 0.35, 
                    label = label_text, 
                    size = 3.0, hjust = 0, vjust = 0,
                    color = "black", fontface = "plain",
                    lineheight = 1.3)
  
  return(p)
}

# ============================================
# 创建三个图
# ============================================

pA <- create_cumulative_plot(
  aj_lc_event,
  "(A) Lung cancer incidence",
  "Cumulative incidence (%)",
  y_max = max(aj_lc_event$curve$cumulative, na.rm = TRUE) * 1.6
)

pB <- create_cumulative_plot(
  aj_lc_death,
  "(B) Lung cancer mortality",
  "Cumulative mortality (%)",
  y_max = max(aj_lc_death$curve$cumulative, na.rm = TRUE) * 1.6
)

pC <- create_cumulative_plot(
  aj_all_death,
  "(C) All-cause mortality",
  "Cumulative mortality (%)",
  y_max = max(aj_all_death$curve$cumulative, na.rm = TRUE) * 1.6
)

# ============================================
# 组合三个图
# ============================================

tryCatch({
  library(cowplot)
  
  # 创建总标题
  title_grob <- ggdraw() + 
    draw_label("Figure 3 Cumulative lung cancer incidence (A), lung cancer mortality (B),\nand all-cause mortality (C) by risk group",
               fontface = "bold", size = 14, hjust = 0.5)
  
  # 组合三个图
  plot_row <- plot_grid(pA, pB, pC, 
                        ncol = 3, 
                        align = "h", axis = "tb")
  
  fig3_combined <- plot_grid(title_grob, plot_row, ncol = 1, rel_heights = c(0.1, 1))
  
  # 导出为TIFF (仅TIFF)
  tiff("/workspace/figure3_three_groups.tiff", width = 20, height = 8, units = "in", res = 300, compression = "lzw")
  print(fig3_combined)
  dev.off()
  
  cat("\nFigure 3 已导出:\n")
  cat("- /workspace/figure3_three_groups.tiff (300 dpi, LZW压缩)\n")
  
}, error = function(e) {
  cat("\ncowplot不可用，使用grid.arrange...\n")
  
  fig3_combined <- grid.arrange(
    pA, pB, pC,
    ncol = 3,
    top = textGrob(
      "Figure 3 Cumulative lung cancer incidence (A), lung cancer mortality (B),\nand all-cause mortality (C) by risk group",
      gp = gpar(fontsize = 14, fontface = "bold")
    )
  )
  
  tiff("/workspace/figure3_three_groups.tiff", width = 20, height = 8, units = "in", res = 300, compression = "lzw")
  grid.draw(fig3_combined)
  dev.off()
  
  cat("\nFigure 3 已导出:\n")
  cat("- /workspace/figure3_three_groups.tiff (300 dpi, LZW压缩)\n")
})

# ============================================
# 导出竞争风险数值表格（CSV）
# ============================================
summary_all <- rbind(
  cbind(Outcome = "Lung cancer incidence", aj_lc_event$summary),
  cbind(Outcome = "Lung cancer mortality", aj_lc_death$summary),
  cbind(Outcome = "All-cause mortality", aj_all_death$summary)
)

fwrite(as.data.table(summary_all), "/workspace/figure3_three_groups_summary.csv")
cat("\n三组竞争风险数值表格已导出:\n")
cat("- /workspace/figure3_three_groups_summary.csv\n")

# ============================================
# 生成宽格式表格（更易于阅读）
# ============================================
summary_wide <- dcast(as.data.table(summary_all), 
                      Outcome + group ~ year, 
                      value.var = "cumulative_rate")
setnames(summary_wide, c("1", "5", "10"), c("Year_1", "Year_5", "Year_10"))

fwrite(summary_wide, "/workspace/figure3_three_groups_summary_wide.csv")
cat("\n宽格式表格已导出:\n")
cat("- /workspace/figure3_three_groups_summary_wide.csv\n")

cat("\n========================================\n")
cat("Figure 3 (三组+竞争风险调整) 导出完成!\n")
cat("========================================\n")
