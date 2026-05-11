# 加载必要的包
library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# 设置中文字体
theme_set(theme_classic(base_size = 12, base_family = "Arial"))

# 从Excel读取数据
library(openxlsx)
df <- read.xlsx("/workspace/表格汇总2026.5.xlsx", sheet = "NCC密度", colNames = FALSE)

# 提取发病率的总体数据 (Total row for Lung cancer incidence)
# 根据数据结构，提取关键行
incidence_row <- 3  # Total for lung cancer incidence
mortality_row <- 27  # Total for lung cancer mortality
allcause_row <- 52  # Total for all-cause mortality

# 解析95% CI的函数
parse_rate_ci <- function(x) {
  if (is.na(x) || x == "-") {
    return(c(NA, NA, NA))
  }
  # 格式: "98.69 (94.59-102.93)"
  x <- as.character(x)
  if (grepl("\\(", x)) {
    rate <- as.numeric(str_extract(x, "^[0-9.]+"))
    ci <- str_extract(x, "\\(([0-9.-]+)-([0-9.-]+)\\)")
    ci_lower <- as.numeric(str_extract(ci, "(?<=\\()[0-9.-]+"))
    ci_upper <- as.numeric(str_extract(ci, "(?<=-)[0-9.-]+(?=\\))"))
    return(c(rate, ci_lower, ci_upper))
  }
  return(c(as.numeric(x), NA, NA))
}

# 提取三组数据
# 列结构: 变量, Overall事件数, Overall人年, Overall率(CI),
#         高危筛查事件数, 高危筛查人年, 高危筛查率(CI),
#         高危未筛查事件数, 高危未筛查人年, 高危未筛查率(CI),
#         低危事件数, 低危人年, 低危率(CI)

# 肺癌发病率 - Total
inc_data <- df[incidence_row, ]
inc_rate_hr_screen <- parse_rate_ci(inc_data[[8]])
inc_rate_hr_unscreen <- parse_rate_ci(inc_data[[12]])
inc_rate_lr <- parse_rate_ci(inc_data[[13]])

# 肺癌死亡率 - Total
mort_data <- df[mortality_row, ]
mort_rate_hr_screen <- parse_rate_ci(mort_data[[8]])
mort_rate_hr_unscreen <- parse_rate_ci(mort_data[[12]])
mort_rate_lr <- parse_rate_ci(mort_data[[13]])

# 全因死亡率 - Total
all_data <- df[allcause_row, ]
all_rate_hr_screen <- parse_rate_ci(all_data[[8]])
all_rate_hr_unscreen <- parse_rate_ci(all_data[[12]])
all_rate_lr <- parse_rate_ci(all_data[[13]])

# 创建数据框
plot_data <- data.frame(
  Outcome = rep(c("Lung Cancer\nIncidence", "Lung Cancer\nMortality", "All-cause\nMortality"), each = 3),
  Group = rep(c("High-risk\nScreened", "High-risk\nUnscreened", "Low-risk\nControl"), times = 3),
  Rate = c(inc_rate_hr_screen[1], inc_rate_hr_unscreen[1], inc_rate_lr[1],
           mort_rate_hr_screen[1], mort_rate_hr_unscreen[1], mort_rate_lr[1],
           all_rate_hr_screen[1], all_rate_hr_unscreen[1], all_rate_lr[1]),
  CI_lower = c(inc_rate_hr_screen[2], inc_rate_hr_unscreen[2], inc_rate_lr[2],
               mort_rate_hr_screen[2], mort_rate_hr_unscreen[2], mort_rate_lr[2],
               all_rate_hr_screen[2], all_rate_hr_unscreen[2], all_rate_lr[2]),
  CI_upper = c(inc_rate_hr_screen[3], inc_rate_hr_unscreen[3], inc_rate_lr[3],
               mort_rate_hr_screen[3], mort_rate_hr_unscreen[3], mort_rate_lr[3],
               all_rate_hr_screen[3], all_rate_hr_unscreen[3], all_rate_lr[3])
)

plot_data$Group <- factor(plot_data$Group, levels = c("High-risk\nScreened", "High-risk\nUnscreened", "Low-risk\nControl"))
plot_data$Outcome <- factor(plot_data$Outcome, levels = c("Lung Cancer\nIncidence", "Lung Cancer\nMortality", "All-cause\nMortality"))

# 定义颜色
colors <- c("#E74C3C", "#3498DB", "#2ECC71")  # 红、蓝、绿

# 创建分组柱状图
p1 <- ggplot(plot_data, aes(x = Outcome, y = Rate, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.8), width = 0.25) +
  scale_fill_manual(values = colors, name = "") +
  labs(y = "Rate per 100,000 person-years", x = "") +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "top",
    legend.box = "horizontal",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# 创建森林图样式的对比图
# 计算相对于低危对照组的RR
rr_data <- plot_data %>%
  group_by(Outcome) %>%
  mutate(
    Ref_rate = Rate[Group == "Low-risk\nControl"],
    RR = Rate / Ref_rate,
    LogRR = log(RR)
  ) %>%
  ungroup()

# 为高危筛查组计算RR
rr_screen <- rr_data[rr_data$Group == "High-risk\nScreened", ]
rr_screen$CI_lower_ratio <- rr_screen$CI_lower / rr_screen$Ref_rate
rr_screen$CI_upper_ratio <- rr_screen$CI_upper / rr_screen$Ref_rate

# 为高危未筛查组计算RR
rr_unscreen <- rr_data[rr_data$Group == "High-risk\nUnscreened", ]
rr_unscreen$CI_lower_ratio <- rr_unscreen$CI_lower / rr_unscreen$Ref_rate
rr_unscreen$CI_upper_ratio <- rr_unscreen$CI_upper / rr_unscreen$Ref_rate

# 合并RR数据
rr_plot_data <- rbind(rr_screen, rr_unscreen)
rr_plot_data <- rr_plot_data[!rr_plot_data$Group == "Low-risk\nControl", ]

# 简化标签
rr_plot_data$Outcome <- gsub("\n", " ", rr_plot_data$Outcome)
rr_plot_data$Group <- gsub("\n", " ", rr_plot_data$Group)

# 创建组合标签
rr_plot_data$Label <- paste0(rr_plot_data$Outcome, "\n(", rr_plot_data$Group, ")")

# 计算RR和CI
rr_plot_data$RR_round <- round(rr_plot_data$RR, 2)
rr_plot_data$CI_text <- paste0("(", round(rr_plot_data$CI_lower_ratio, 2), "-", round(rr_plot_data$CI_upper_ratio, 2), ")")

# 创建森林图
p2 <- ggplot(rr_plot_data, aes(x = Label, y = RR, color = Group)) +
  geom_point(position = position_dodge(width = 0.5), size = 4, shape = 18) +
  geom_errorbar(aes(ymin = CI_lower_ratio, ymax = CI_upper_ratio),
                position = position_dodge(width = 0.5), width = 0.15, size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", size = 1) +
  scale_color_manual(values = c("#E74C3C", "#3498DB"), name = "") +
  labs(y = "Rate Ratio (vs Low-risk Control)", x = "", title = "Rate Ratio with 95% CI") +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# 添加RR数值标注
p2 <- p2 + geom_text(aes(y = RR + 0.15, label = paste0(RR_round, "\n", CI_text)),
                     position = position_dodge(width = 0.5), size = 3, vjust = 0)

# 组合两个图
combined_plot <- ggarrange(
  p1 + ggtitle("A. Event Rates by Group") + theme(plot.title = element_text(face = "bold", hjust = 0)),
  p2 + theme(legend.position = "none") + ggtitle("B. Rate Ratio (vs Low-risk Control)"),
  ncol = 2, widths = c(1.3, 1), common.legend = TRUE, legend = "top"
)

# 添加总标题
combined_plot <- annotate_figure(combined_plot, top = text_grob("Comparison of Lung Cancer Outcomes Across Risk Groups", size = 14, face = "bold"))

# 保存图片
ggsave("/workspace/figure_comparison_rates.png", combined_plot, width = 14, height = 7, dpi = 300, bg = "white")

# 打印数据汇总
cat("\n===== 数据汇总 =====\n")
print(plot_data)
cat("\n===== RR (Rate Ratio) vs Low-risk Control =====\n")
print(rr_plot_data[, c("Outcome", "Group", "Rate", "Ref_rate", "RR", "CI_lower_ratio", "CI_upper_ratio")])

cat("\n图片已保存至: /workspace/figure_comparison_rates.png\n")
