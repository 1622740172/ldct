library(openxlsx)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(stringr)

theme_set(theme_classic(base_size = 12))

df <- read.xlsx("/workspace/表格汇总2026.5.xlsx", sheet = "NCC密度", colNames = FALSE)

parse_rate_ci <- function(x) {
  if (is.na(x) || x == "-" || is.null(x) || is.nan(x)) return(c(NA, NA, NA))
  x <- as.character(x)
  if (!grepl("(", x, fixed = TRUE)) {
    rate <- as.numeric(x)
    if (is.na(rate)) return(c(NA, NA, NA))
    return(c(rate, NA, NA))
  }
  rate <- as.numeric(str_extract(x, "^[0-9]+\\.[0-9]+"))
  ci_parts <- str_match(x, "\\(([0-9]+\\.[0-9]+)-([0-9]+\\.[0-9]+)\\)")
  ci_lower <- as.numeric(ci_parts[2])
  ci_upper <- as.numeric(ci_parts[3])
  return(c(rate, ci_lower, ci_upper))
}

cat("=== NCC密度 数据提取 ===\n")

inc_row <- 4
mort_row <- 28
all_row <- 52

inc_hr_screen <- parse_rate_ci(df[inc_row, 7])
inc_hr_unscreen <- parse_rate_ci(df[inc_row, 10])
inc_lr <- parse_rate_ci(df[inc_row, 13])

mort_hr_screen <- parse_rate_ci(df[mort_row, 7])
mort_hr_unscreen <- parse_rate_ci(df[mort_row, 10])
mort_lr <- parse_rate_ci(df[mort_row, 13])

all_hr_screen <- parse_rate_ci(df[all_row, 7])
all_hr_unscreen <- parse_rate_ci(df[all_row, 10])
all_lr <- parse_rate_ci(df[all_row, 13])

cat("\n=== 肺癌发病率 (每10万人年) ===\n")
cat("高危筛查组: ", inc_hr_screen[1], " (", inc_hr_screen[2], "-", inc_hr_screen[3], ")\n", sep="")
cat("高危未筛查组: ", inc_hr_unscreen[1], " (", inc_hr_unscreen[2], "-", inc_hr_unscreen[3], ")\n", sep="")
cat("低危对照组: ", inc_lr[1], " (", inc_lr[2], "-", inc_lr[3], ")\n", sep="")

cat("\n=== 肺癌死亡率 (每10万人年) ===\n")
cat("高危筛查组: ", mort_hr_screen[1], " (", mort_hr_screen[2], "-", mort_hr_screen[3], ")\n", sep="")
cat("高危未筛查组: ", mort_hr_unscreen[1], " (", mort_hr_unscreen[2], "-", mort_hr_unscreen[3], ")\n", sep="")
cat("低危对照组: ", mort_lr[1], " (", mort_lr[2], "-", mort_lr[3], ")\n", sep="")

cat("\n=== 全因死亡率 (每10万人年) ===\n")
cat("高危筛查组: ", all_hr_screen[1], " (", all_hr_screen[2], "-", all_hr_screen[3], ")\n", sep="")
cat("高危未筛查组: ", all_hr_unscreen[1], " (", all_hr_unscreen[2], "-", all_hr_unscreen[3], ")\n", sep="")
cat("低危对照组: ", all_lr[1], " (", all_lr[2], "-", all_lr[3], ")\n", sep="")

plot_data <- data.frame(
  Outcome = rep(c("Lung Cancer\nIncidence", "Lung Cancer\nMortality", "All-cause\nMortality"), each = 3),
  Group = rep(c("High-risk\nScreened", "High-risk\nUnscreened", "Low-risk\nControl"), times = 3),
  Rate = c(inc_hr_screen[1], inc_hr_unscreen[1], inc_lr[1],
           mort_hr_screen[1], mort_hr_unscreen[1], mort_lr[1],
           all_hr_screen[1], all_hr_unscreen[1], all_lr[1]),
  CI_lower = c(inc_hr_screen[2], inc_hr_unscreen[2], inc_lr[2],
               mort_hr_screen[2], mort_hr_unscreen[2], mort_lr[2],
               all_hr_screen[2], all_hr_unscreen[2], all_lr[2]),
  CI_upper = c(inc_hr_screen[3], inc_hr_unscreen[3], inc_lr[3],
               mort_hr_screen[3], mort_hr_unscreen[3], mort_lr[3],
               all_hr_screen[3], all_hr_unscreen[3], all_lr[3])
)

plot_data$Group <- factor(plot_data$Group, levels = c("High-risk\nScreened", "High-risk\nUnscreened", "Low-risk\nControl"))
plot_data$Outcome <- factor(plot_data$Outcome, levels = c("Lung Cancer\nIncidence", "Lung Cancer\nMortality", "All-cause\nMortality"))

colors <- c("#E74C3C", "#3498DB", "#2ECC71")

p1 <- ggplot(plot_data, aes(x = Outcome, y = Rate, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.8) +
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

rr_data <- plot_data %>%
  group_by(Outcome) %>%
  mutate(
    Ref_rate = Rate[Group == "Low-risk\nControl"],
    RR = Rate / Ref_rate,
    RR_CI_lower = CI_lower / Ref_rate,
    RR_CI_upper = CI_upper / Ref_rate
  ) %>%
  ungroup()

rr_plot_data <- rr_data[!rr_data$Group == "Low-risk\nControl", ]
rr_plot_data$Label <- paste0(rr_plot_data$Outcome, " (", rr_plot_data$Group, ")")
rr_plot_data$Label <- gsub("\n", " ", rr_plot_data$Label)
rr_plot_data$Group <- gsub("\n", " ", rr_plot_data$Group)

y_max <- max(rr_plot_data$RR_CI_upper, na.rm = TRUE) * 1.3

p2 <- ggplot(rr_plot_data, aes(x = Label, y = RR, color = Group, shape = Group)) +
  geom_point(position = position_dodge(width = 0.5), size = 5) +
  geom_errorbar(aes(ymin = RR_CI_lower, ymax = RR_CI_upper),
                position = position_dodge(width = 0.5), width = 0.15, linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 1) +
  scale_color_manual(values = c("#E74C3C", "#3498DB"), name = "") +
  scale_shape_manual(values = c(18, 16), name = "") +
  labs(y = "Rate Ratio (vs Low-risk Control)", x = "", title = "B. Rate Ratio (95% CI)") +
  ylim(0, y_max) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 15, hjust = 1),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0)
  )

combined_plot <- ggarrange(
  p1 + ggtitle("A. Event Rates by Group") + theme(plot.title = element_text(face = "bold", hjust = 0)),
  p2,
  ncol = 2, widths = c(1.2, 1.1), common.legend = TRUE, legend = "top"
)

combined_plot <- annotate_figure(combined_plot, top = text_grob("Comparison of Outcomes: High-risk Screened vs Unscreened vs Low-risk Control", size = 14, face = "bold"))

ggsave("/workspace/figure_comparison_rates.png", combined_plot, width = 14, height = 7, dpi = 300, bg = "white")

cat("\n===== Rate Ratio (RR) vs Low-risk Control =====\n")
print(rr_plot_data[, c("Outcome", "Group", "Rate", "Ref_rate", "RR", "RR_CI_lower", "RR_CI_upper")])

cat("\n图片已保存至: /workspace/figure_comparison_rates.png\n")
