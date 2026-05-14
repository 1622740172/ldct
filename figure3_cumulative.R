# ============================================
# Figure 3: Cumulative incidence/mortality curves
# (A) Lung cancer incidence
# (B) Lung cancer mortality  
# (C) All-cause mortality
# Comparing screened (as_lc_risk=1) vs non-screened (as_lc_risk=2)
# X-axis: years, 1-year intervals
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
# 计算累积发生率 (Aalen-Johansen / 竞争风险方法)
# 使用生存分析中的累积发生率函数
# ============================================

calculate_cumulative <- function(dt, event_col, time_col, max_years = 12) {
  # 按组分别计算
  groups <- unique(dt$group)
  result <- data.frame()
  
  for (g in groups) {
    sub_dt <- dt[group == g]
    n_total <- nrow(sub_dt)
    
    # 创建时间网格 (0到max_years, 步长0.1年用于平滑曲线)
    time_grid <- seq(0, max_years, by = 0.1)
    
    cum_inc <- numeric(length(time_grid))
    
    for (i in 1:length(time_grid)) {
      t <- time_grid[i]
      # 在时间点t之前发生事件的人数
      events_by_t <- sub_dt[get(time_col) <= t & get(event_col) == 1, .N]
      # 在时间点t仍在随访中的人数 (包括t之后失访和t之后发生事件的)
      at_risk <- sub_dt[get(time_col) >= t, .N]
      
      if (n_total > 0) {
        cum_inc[i] <- events_by_t / n_total * 100  # 转换为百分比
      } else {
        cum_inc[i] <- 0
      }
    }
    
    result <- rbind(result, data.frame(
      time = time_grid,
      cumulative = cum_inc,
      group = g,
      stringsAsFactors = FALSE
    ))
  }
  
  return(result)
}

# 更精确的方法：使用生存包计算累积发生率
# 安装survival包
tryCatch({
  library(survival)
}, error = function(e) {
  install.packages("survival", repos = "https://cloud.r-project.org/")
  library(survival)
})

# ============================================
# 使用survfit计算生存曲线，然后转换为累积发生率
# ============================================

calculate_cumulative_surv <- function(dt, event_col, time_col, max_years = 12) {
  groups <- unique(dt$group)
  result <- data.frame()
  
  for (g in groups) {
    sub_dt <- dt[group == g]
    
    # 创建生存对象
    # status: 1=event, 0=censored
    surv_obj <- Surv(time = sub_dt[[time_col]], event = sub_dt[[event_col]])
    
    # 拟合Kaplan-Meier生存曲线
    fit <- survfit(surv_obj ~ 1, data = sub_dt)
    
    # 提取时间点和生存概率
    time_points <- fit$time / 365.25  # 转换为年
    surv_prob <- fit$surv
    
    # 计算累积发生率 = 1 - 生存概率
    cum_inc <- (1 - surv_prob) * 100  # 转换为百分比
    
    # 限制在max_years以内
    valid_idx <- time_points <= max_years
    time_points <- time_points[valid_idx]
    cum_inc <- cum_inc[valid_idx]
    
    # 添加起点 (0, 0)
    time_points <- c(0, time_points)
    cum_inc <- c(0, cum_inc)
    
    # 插值到更细的时间网格用于平滑曲线
    time_grid <- seq(0, max_years, by = 0.05)
    cum_inc_interp <- approx(time_points, cum_inc, xout = time_grid, method = "constant", f = 0)$y
    cum_inc_interp[is.na(cum_inc_interp)] <- max(cum_inc_interp, na.rm = TRUE)
    
    result <- rbind(result, data.frame(
      time = time_grid,
      cumulative = cum_inc_interp,
      group = g,
      stringsAsFactors = FALSE
    ))
  }
  
  return(result)
}

# ============================================
# 计算三个结局的累积曲线
# ============================================

cat("计算肺癌发病率累积曲线...\n")
cum_lc_event <- calculate_cumulative_surv(df_high, "LC_event", "cancer_time", max_years = 12)

cat("计算肺癌死亡率累积曲线...\n")
cum_lc_death <- calculate_cumulative_surv(df_high, "LC_death", "death_time", max_years = 12)

cat("计算全因死亡率累积曲线...\n")
cum_all_death <- calculate_cumulative_surv(df_high, "all_cause_death", "death_time", max_years = 12)

# ============================================
# 设置颜色和样式
# ============================================
group_colors <- c("Screened" = "#2E86AB", "Non-screened" = "#A23B72")
group_linetypes <- c("Screened" = "solid", "Non-screened" = "dashed")

# ============================================
# 创建累积曲线图函数
# ============================================

create_cumulative_plot <- function(cum_data, title_text, y_label, y_max = NULL) {
  if (is.null(y_max)) {
    y_max <- max(cum_data$cumulative, na.rm = TRUE) * 1.2
  }
  
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
      legend.position = c(0.25, 0.85),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.background = element_rect(fill = "white", color = "gray80"),
      legend.box.background = element_rect(color = "gray80"),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  return(p)
}

# ============================================
# 创建三个图
# ============================================

pA <- create_cumulative_plot(
  cum_lc_event,
  "(A) Lung cancer incidence",
  "Cumulative incidence (%)",
  y_max = max(cum_lc_event$cumulative, na.rm = TRUE) * 1.5
)

pB <- create_cumulative_plot(
  cum_lc_death,
  "(B) Lung cancer mortality",
  "Cumulative mortality (%)",
  y_max = max(cum_lc_death$cumulative, na.rm = TRUE) * 1.5
)

pC <- create_cumulative_plot(
  cum_all_death,
  "(C) All-cause mortality",
  "Cumulative mortality (%)",
  y_max = max(cum_all_death$cumulative, na.rm = TRUE) * 1.5
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
  tiff("/workspace/figure3_cumulative.tiff", width = 18, height = 7, units = "in", res = 300, compression = "lzw")
  print(fig3_combined)
  dev.off()
  
  cat("\nFigure 3 已导出:\n")
  cat("- /workspace/figure3_cumulative.tiff (300 dpi, LZW压缩)\n")
  
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
  
  tiff("/workspace/figure3_cumulative.tiff", width = 18, height = 7, units = "in", res = 300, compression = "lzw")
  grid.draw(fig3_combined)
  dev.off()
  
  cat("\nFigure 3 已导出:\n")
  cat("- /workspace/figure3_cumulative.tiff (300 dpi, LZW压缩)\n")
})

cat("\n========================================\n")
cat("Figure 3 导出完成!\n")
cat("========================================\n")
