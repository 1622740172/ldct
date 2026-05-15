# ============================================
# Figure 3: Cumulative incidence/mortality curves with competing risk adjustment
# (A) Lung cancer incidence
# (B) Lung cancer mortality  
# (C) All-cause mortality
# Comparing screened (as_lc_risk=1) vs non-screened (as_lc_risk=2)
# X-axis: years, 1-year intervals
# Using Aalen-Johansen estimator for competing risks
# ============================================

library(data.table)
library(ggplot2)
library(cowplot)

# 读取数据
df <- fread("/workspace/匹配结果_最终修正版.csv", encoding="UTF-8")

# 筛选高危人群 (screened=1 vs non-screened=2)
df_high <- df[as_lc_risk %in% c(1, 2)]
df_high[, group := ifelse(as_lc_risk == 1, "Screened", "Non-screened")]

# 转换时间为年
df_high[, time_years := cancer_time / 365.25]

# 确保事件变量是数值型
df_high[, LC_event := as.numeric(LC_event)]
df_high[, LC_death := as.numeric(LC_death)]
df_high[, all_cause_death := as.numeric(all_cause_death)]

# ============================================
# 竞争风险模型：Aalen-Johansen 估计器
# 对于肺癌发病：竞争事件是死亡（先死亡则无法观察到肺癌发病）
# 对于肺癌死亡：竞争事件是非肺癌死亡
# 对于全因死亡：无竞争事件（使用KM即可）
# ============================================

# 手动实现Aalen-Johansen估计器
calculate_aj_cumulative <- function(dt, event_col, competing_col, time_col, max_years = 12) {
  groups <- unique(dt$group)
  result <- data.frame()
  summary_table <- data.frame()
  
  for (g in groups) {
    sub_dt <- dt[group == g]
    n_total <- nrow(sub_dt)
    
    # 获取所有唯一的事件时间（按天）
    all_times <- sort(unique(sub_dt[[time_col]]))
    all_times <- all_times[all_times > 0]
    
    # 初始化
    cum_inc <- 0
    surv_overall <- 1
    
    # 存储结果
    time_points <- c(0)
    cum_values <- c(0)
    
    for (t in all_times) {
      # 在时间点t的风险集人数
      at_risk <- sum(sub_dt[[time_col]] >= t)
      
      if (at_risk == 0) next
      
      # 在时间点t发生目标事件的人数
      d_event <- sum(sub_dt[[time_col]] == t & sub_dt[[event_col]] == 1)
      
      # 在时间点t发生竞争事件的人数
      d_competing <- sum(sub_dt[[time_col]] == t & sub_dt[[competing_col]] == 1)
      
      # 总事件数
      d_total <- d_event + d_competing
      
      if (at_risk > 0) {
        # 更新总体生存概率
        surv_overall <- surv_overall * (1 - d_total / at_risk)
        
        # 更新累积发生率
        if (at_risk > 0) {
          cum_inc <- cum_inc + surv_overall * (d_event / at_risk)
        }
      }
      
      time_points <- c(time_points, t / 365.25)
      cum_values <- c(cum_values, cum_inc * 100)  # 转换为百分比
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
# 计算三个结局的累积曲线（竞争风险调整）
# ============================================

cat("计算肺癌发病率累积曲线（竞争风险调整）...\n")
# 对于肺癌发病，竞争事件是：先发生死亡（LC_death或non_lc_death）
# 使用 all_cause_death 作为竞争事件（如果先死亡则无法观察到肺癌发病）
aj_lc_event <- calculate_aj_cumulative(
  df_high, 
  event_col = "LC_event", 
  competing_col = "all_cause_death",  # 死亡是竞争事件
  time_col = "cancer_time",
  max_years = 12
)

cat("计算肺癌死亡率累积曲线（竞争风险调整）...\n")
# 对于肺癌死亡，竞争事件是：非肺癌死亡
# 创建non_lc_death变量
df_high[, non_lc_death := ifelse(all_cause_death == 1 & LC_death == 0, 1, 0)]
aj_lc_death <- calculate_aj_cumulative(
  df_high, 
  event_col = "LC_death", 
  competing_col = "non_lc_death",  # 非肺癌死亡是竞争事件
  time_col = "death_time",
  max_years = 12
)

cat("计算全因死亡率累积曲线...\n")
# 全因死亡无竞争事件，使用KM方法
aj_all_death <- calculate_aj_cumulative(
  df_high, 
  event_col = "all_cause_death", 
  competing_col = "LC_event",  # 占位，实际上不会用到
  time_col = "death_time",
  max_years = 12
)

# ============================================
# 设置颜色和样式
# ============================================
group_colors <- c("Screened" = "#2E86AB", "Non-screened" = "#A23B72")
group_linetypes <- c("Screened" = "solid", "Non-screened" = "dashed")

# ============================================
# 创建累积曲线图函数（带数值标注）
# ============================================

create_cumulative_plot <- function(aj_result, title_text, y_label, y_max = NULL) {
  cum_data <- aj_result$curve
  summary_data <- aj_result$summary
  
  if (is.null(y_max)) {
    y_max <- max(cum_data$cumulative, na.rm = TRUE) * 1.5
  }
  
  # 准备标注文本
  label_text <- ""
  for (g in unique(summary_data$group)) {
    sub <- summary_data[summary_data$group == g, ]
    label_text <- paste0(label_text, g, ": ")
    rates <- paste0(sub$year, "y=", sub$cumulative_rate, "%")
    label_text <- paste0(label_text, paste(rates, collapse = ", "), "\n")
  }
  label_text <- trimws(label_text)
  
  p <- ggplot(cum_data, aes(x = time, y = cumulative, color = group, linetype = group)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = group_colors) +
    scale_linetype_manual(values = group_linetypes) +
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
      legend.position = "none",  # 隐藏默认图例，使用自定义标注
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  # 添加数值标注在图内部
  p <- p + annotate("text", x = 7, y = y_max * 0.9, 
                    label = label_text, 
                    size = 3.2, hjust = 0, vjust = 1,
                    color = "black", fontface = "plain")
  
  return(p)
}

# ============================================
# 创建三个图
# ============================================

pA <- create_cumulative_plot(
  aj_lc_event,
  "(A) Lung cancer incidence",
  "Cumulative incidence (%)",
  y_max = max(aj_lc_event$curve$cumulative, na.rm = TRUE) * 1.8
)

pB <- create_cumulative_plot(
  aj_lc_death,
  "(B) Lung cancer mortality",
  "Cumulative mortality (%)",
  y_max = max(aj_lc_death$curve$cumulative, na.rm = TRUE) * 1.8
)

pC <- create_cumulative_plot(
  aj_all_death,
  "(C) All-cause mortality",
  "Cumulative mortality (%)",
  y_max = max(aj_all_death$curve$cumulative, na.rm = TRUE) * 1.8
)

# ============================================
# 组合三个图
# ============================================

tryCatch({
  library(cowplot)
  
  # 创建总标题
  title_grob <- ggdraw() + 
    draw_label("Figure 3 Cumulative lung cancer incidence (A), lung cancer mortality (B),\nand all-cause mortality (C)",
               fontface = "bold", size = 14, hjust = 0.5)
  
  # 组合三个图
  plot_row <- plot_grid(pA, pB, pC, 
                        ncol = 3, 
                        align = "h", axis = "tb")
  
  fig3_combined <- plot_grid(title_grob, plot_row, ncol = 1, rel_heights = c(0.1, 1))
  
  # 导出为TIFF (仅TIFF)
  tiff("/workspace/figure3_cumulative_cr.tiff", width = 18, height = 7, units = "in", res = 300, compression = "lzw")
  print(fig3_combined)
  dev.off()
  
  cat("\nFigure 3 已导出:\n")
  cat("- /workspace/figure3_cumulative_cr.tiff (300 dpi, LZW压缩)\n")
  
}, error = function(e) {
  cat("\ncowplot不可用，使用grid.arrange...\n")
  
  fig3_combined <- grid.arrange(
    pA, pB, pC,
    ncol = 3,
    top = textGrob(
      "Figure 3 Cumulative lung cancer incidence (A), lung cancer mortality (B),\nand all-cause mortality (C)",
      gp = gpar(fontsize = 14, fontface = "bold")
    )
  )
  
  tiff("/workspace/figure3_cumulative_cr.tiff", width = 18, height = 7, units = "in", res = 300, compression = "lzw")
  grid.draw(fig3_combined)
  dev.off()
  
  cat("\nFigure 3 已导出:\n")
  cat("- /workspace/figure3_cumulative_cr.tiff (300 dpi, LZW压缩)\n")
})

# ============================================
# 导出竞争风险数值表格（CSV）
# ============================================
summary_all <- rbind(
  cbind(Outcome = "Lung cancer incidence", aj_lc_event$summary),
  cbind(Outcome = "Lung cancer mortality", aj_lc_death$summary),
  cbind(Outcome = "All-cause mortality", aj_all_death$summary)
)

fwrite(as.data.table(summary_all), "/workspace/figure3_competing_risk_summary.csv")
cat("\n竞争风险数值表格已导出:\n")
cat("- /workspace/figure3_competing_risk_summary.csv\n")

cat("\n========================================\n")
cat("Figure 3 (竞争风险调整) 导出完成!\n")
cat("========================================\n")
