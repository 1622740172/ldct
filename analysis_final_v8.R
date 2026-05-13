# ============================================
# 最终版v8：修正人数和RR显示
# 左侧表格：Subgroup + 三个结局的RR(95%CI) + P值
# 右侧：三个纯森林图
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
# 计算RR和详细统计量
# ============================================
calculate_rr <- function(dt, outcome_col, time_col) {
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
    n_screened = n1, n_unscreened = n0,
    events_screened = e1, events_unscreened = e0,
    py_screened = py1, py_unscreened = py0,
    rate_screened = rate1, rate_unscreened = rate0,
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
        n_screened = NA, n_unscreened = NA,
        events_screened = NA, events_unscreened = NA,
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
      
      cat(sprintf("%s: n_screened=%d, n_unscreened=%d, RR=%s %s, P=%s\n", 
                  label, res$n_screened, res$n_unscreened, rr_str, ci_str, p_str))
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
  
  rr_col <- paste0("rr_", outcome)
  ci_lower_col <- paste0("ci_lower_", outcome)
  ci_upper_col <- paste0("ci_upper_", outcome)
  p_col <- paste0("p_", outcome)
  n_s_col <- paste0("n_s_", outcome)
  n_u_col <- paste0("n_u_", outcome)
  
  plot_df[[rr_col]] <- NA
  plot_df[[ci_lower_col]] <- NA
  plot_df[[ci_upper_col]] <- NA
  plot_df[[p_col]] <- NA
  plot_df[[n_s_col]] <- NA
  plot_df[[n_u_col]] <- NA
  
  for (i in 1:n) {
    res <- res_list[[i]]
    if (res$type == "data") {
      plot_df[[rr_col]][i] <- res$rr
      plot_df[[ci_lower_col]][i] <- res$ci_lower
      plot_df[[ci_upper_col]][i] <- res$ci_upper
      plot_df[[p_col]][i] <- res$p_value
      plot_df[[n_s_col]][i] <- res$n_screened
      plot_df[[n_u_col]][i] <- res$n_unscreened
    }
  }
}

# ============================================
# 创建左侧表格Grob - 显示三个结局的RR
# ============================================

create_left_table <- function(plot_df) {
  n <- nrow(plot_df)
  
  # 创建表格数据 - 每个结局一列RR
  table_data <- data.frame(
    Group = character(n),
    RR_Incidence = character(n),
    RR_Mortality = character(n),
    RR_AllCause = character(n),
    P_Incidence = character(n),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n) {
    if (plot_df$type[i] == "header") {
      table_data$Group[i] <- plot_df$display_label[i]
      table_data$RR_Incidence[i] <- ""
      table_data$RR_Mortality[i] <- ""
      table_data$RR_AllCause[i] <- ""
      table_data$P_Incidence[i] <- ""
    } else {
      table_data$Group[i] <- plot_df$display_label[i]
      
      # 肺癌发病率 RR
      rr <- plot_df$rr_LC_event[i]
      ci_l <- plot_df$ci_lower_LC_event[i]
      ci_u <- plot_df$ci_upper_LC_event[i]
      if (!is.na(rr) && !is.infinite(rr)) {
        table_data$RR_Incidence[i] <- sprintf("%.2f (%.2f-%.2f)", rr, ci_l, ci_u)
      } else {
        table_data$RR_Incidence[i] <- "N/A"
      }
      
      # 肺癌死亡率 RR
      rr <- plot_df$rr_LC_death[i]
      ci_l <- plot_df$ci_lower_LC_death[i]
      ci_u <- plot_df$ci_upper_LC_death[i]
      if (!is.na(rr) && !is.infinite(rr)) {
        table_data$RR_Mortality[i] <- sprintf("%.2f (%.2f-%.2f)", rr, ci_l, ci_u)
      } else {
        table_data$RR_Mortality[i] <- "N/A"
      }
      
      # 全因死亡率 RR
      rr <- plot_df$rr_all_cause_death[i]
      ci_l <- plot_df$ci_lower_all_cause_death[i]
      ci_u <- plot_df$ci_upper_all_cause_death[i]
      if (!is.na(rr) && !is.infinite(rr)) {
        table_data$RR_AllCause[i] <- sprintf("%.2f (%.2f-%.2f)", rr, ci_l, ci_u)
      } else {
        table_data$RR_AllCause[i] <- "N/A"
      }
      
      # P值（使用肺癌发病率的P值）
      p <- plot_df$p_LC_event[i]
      if (!is.na(p)) {
        if (p < 0.0001) {
          table_data$P_Incidence[i] <- "<0.0001"
        } else {
          table_data$P_Incidence[i] <- sprintf("%.4f", p)
        }
      } else {
        table_data$P_Incidence[i] <- "N/A"
      }
    }
  }
  
  # 创建tableGrob
  tg <- tableGrob(
    d = table_data,
    rows = NULL,
    cols = c("Subgroup", "LC incidence\nRR (95% CI)", "LC mortality\nRR (95% CI)", "All-cause mortality\nRR (95% CI)", "P"),
    theme = ttheme_minimal(
      core = list(
        fg_params = list(fontsize = 7.5, fontface = "plain", col = "black", lineheight = 0.85),
        bg_params = list(fill = "white", col = NA),
        padding = unit(c(1.5, 2), "mm")
      ),
      colhead = list(
        fg_params = list(fontsize = 8, fontface = "bold", col = "black"),
        bg_params = list(fill = "gray90", col = "gray70"),
        padding = unit(c(2, 2), "mm")
      )
    )
  )
  
  # 设置列宽
  tg$widths <- unit(c(0.18, 0.22, 0.22, 0.22, 0.16), "npc")
  
  # 为标题行添加粗体
  for (i in 1:n) {
    if (plot_df$type[i] == "header") {
      for (j in 1:5) {
        idx <- (i) * 5 + j
        if (idx <= length(tg$grobs)) {
          tg$grobs[[idx]]$gp <- gpar(fontface = "bold", fontsize = 8, col = "black")
        }
      }
    }
  }
  
  return(tg)
}

# ============================================
# 创建右侧森林图
# ============================================

create_forest_plot <- function(plot_df, outcome_name) {
  color <- outcome_colors[outcome_name]
  
  rr_col <- paste0("rr_", outcome_name)
  ci_lower_col <- paste0("ci_lower_", outcome_name)
  ci_upper_col <- paste0("ci_upper_", outcome_name)
  
  valid_idx <- !is.na(plot_df[[rr_col]]) & !is.infinite(plot_df[[rr_col]])
  plot_df_valid <- plot_df[valid_idx, ]
  
  all_ci <- c(plot_df[[ci_lower_col]], plot_df[[ci_upper_col]])
  all_ci <- all_ci[!is.na(all_ci) & !is.infinite(all_ci)]
  x_min <- min(0.2, min(all_ci) * 0.8)
  x_max <- max(5, max(all_ci) * 1.2)
  
  p <- ggplot() +
    geom_segment(data = plot_df_valid,
                 aes(x = .data[[ci_lower_col]], xend = .data[[ci_upper_col]], y = y_pos, yend = y_pos),
                 color = color, linewidth = 1.5, alpha = 0.7) +
    geom_point(data = plot_df_valid,
               aes(x = .data[[rr_col]], y = y_pos),
               color = color, size = 4, shape = 21, fill = color, stroke = 0.8) +
    geom_vline(xintercept = 1, color = "black", linewidth = 0.8) +
    scale_x_log10(limits = c(x_min, x_max),
                  breaks = c(0.5, 1, 1.5, 2, 2.5),
                  labels = c("0.5", "1", "1.5", "2", "2.5")) +
    scale_y_continuous(breaks = plot_df$y_pos, 
                       labels = NULL,
                       expand = expansion(add = c(0.5, 0.5))) +
    labs(title = outcome_labels[outcome_name], x = "RR (95% CI)", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(color = color, size = 13, face = "bold", hjust = 0.5),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p)
}

# ============================================
# 创建所有组件并组合
# ============================================

tg <- create_left_table(plot_df)
p1 <- create_forest_plot(plot_df, "LC_event")
p2 <- create_forest_plot(plot_df, "LC_death")
p3 <- create_forest_plot(plot_df, "all_cause_death")

# 组合
final_plot <- grid.arrange(
  tg,
  arrangeGrob(
    p1 + theme(plot.margin = margin(10, 5, 10, 5)),
    p2 + theme(plot.margin = margin(10, 5, 10, 5)),
    p3 + theme(plot.margin = margin(10, 5, 10, 5)),
    ncol = 3
  ),
  ncol = 2,
  widths = c(0.42, 0.58),
  top = textGrob(
    "Figure 1 The crude rate ratios of lung cancer incidence density, lung cancer mortality,\nand all-cause mortality comparing the screened with non-screened groups",
    gp = gpar(fontsize = 14, fontface = "bold")
  )
)

# 保存为TIFF
tiff("/workspace/forest_plot_final_v8.tiff", width = 32, height = 16, units = "in", res = 300, compression = "lzw")
grid.draw(final_plot)
dev.off()

png("/workspace/forest_plot_final_v8.png", width = 3200, height = 1600, res = 150)
grid.draw(final_plot)
dev.off()

cat("\n图片已保存:\n")
cat("- /workspace/forest_plot_final_v8.tiff (300 dpi, LZW压缩)\n")
cat("- /workspace/forest_plot_final_v8.png (预览)\n")

# 导出结果表格
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
        RR = ifelse(is.na(res$rr) || is.infinite(res$rr), NA, round(res$rr, 2)),
        CI_Lower = ifelse(is.na(res$ci_lower), NA, round(res$ci_lower, 2)),
        CI_Upper = ifelse(is.na(res$ci_upper), NA, round(res$ci_upper, 2)),
        P_Value = ifelse(is.na(res$p_value), NA, res$p_value)
      ))
    }
  }
}

fwrite(result_table, "/workspace/rr_results_final_v8.csv")
cat("结果表格已保存到 /workspace/rr_results_final_v8.csv\n")
