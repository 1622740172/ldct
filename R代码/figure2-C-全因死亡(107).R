# ===================== 0. 环境彻底净化 =====================
rm(list = ls())
gc()


# ==============================================================================
# 修正版代码③：全因死亡 K-M 曲线 + 绝对准确的风险表 (增加10年截断)
# ==============================================================================

# 1. 加载必需包与环境设置
library(survival)
library(dplyr)
library(ggplot2)
library(data.table)
library(showtext)
library(patchwork) 

# 字体设置（确保中文显示）
showtext_auto()
font_add("SimHei", "simhei.ttf") 
showtext_opts(dpi = 300) 

# 2. 数据读取与列名配置
input_file_path <- "D:/匹配结果_最终修正版.csv"

if (!file.exists(input_file_path)) {
  stop("找不到文件，请检查 D 盘路径是否正确或文件名是否有误。")
}

cat("正在读取数据...\n")
raw_data <- fread(input_file_path, header = TRUE, quote = "", fill = TRUE, 
                  data.table = FALSE, na.strings = c("", "NA", "NULL"))

# 3. 数据清洗与状态定义 (⭐⭐ 核心修改：加入 10 年截断逻辑 ⭐⭐)
data_all_cause <- raw_data %>%
  select(
    death_time = death_time,         # 全因死亡统一看 death_time
    all_death  = all_cause_death,    # 全因死亡事件列
    group      = as_lc_risk
  ) %>%
  mutate(
    # 直接将天数除以 365.25 转换为年
    time_years = as.numeric(as.character(death_time)) / 365.25,
    status     = as.numeric(as.character(all_death)),   # 1=死亡，0=存活/删失
    group      = as.numeric(as.character(group))
  ) %>%
  # -------------------- 新增：十年截断处理 --------------------
# 如果随访时间大于10年，则第10年时状态转为删失(0)，随访时间截断为10
mutate(
  status     = ifelse(time_years > 10, 0, status),
  time_years = ifelse(time_years > 10, 10, time_years)
) %>%
  # ------------------------------------------------------------
filter(!is.na(time_years) & time_years > 0 & !is.na(group))

# 定义分组标签
data_all_cause$risk_group <- factor(
  data_all_cause$group,
  levels = c(0, 1, 2),
  labels = c("低危未筛查组", "高危筛查组", "高危未筛查组") 
)

# 4. 统计分析：K-M 拟合与 Log-rank 检验
fit_km <- survfit(Surv(time_years, status) ~ risk_group, data = data_all_cause)

# 计算 Log-rank P值
log_rank_test <- survdiff(Surv(time_years, status) ~ risk_group, data = data_all_cause)
p_value <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)
p_text <- ifelse(p_value < 0.001, 
                 "Log-rank 检验 P < 0.001", 
                 paste0("Log-rank 检验 P = ", round(p_value, 3)))

# 5. 提取绘图数据（累积死亡率）
plot_data <- data.frame(
  time     = fit_km$time,
  surv     = fit_km$surv,
  group    = rep(names(fit_km$strata), fit_km$strata)
) %>%
  mutate(
    cum_mortality = (1 - surv) * 100000,   # 转换为 1/100,000
    group = gsub("risk_group=", "", group)
  )

plot_data$group <- factor(plot_data$group, 
                          levels = c("低危未筛查组", "高危筛查组", "高危未筛查组"))

# 6. 绘图配置 (⭐⭐ 核心修改：时间刻度和最大X轴调整为 10 ⭐⭐)
group_colors <- c("低危未筛查组" = "#95e1d3", 
                  "高危筛查组"   = "#f38181", 
                  "高危未筛查组" = "#fce38a")

# 取消 12.5，将最大展示刻度固定为 10
target_times <- c(0, 2.5, 5, 7.5, 10)
max_x <- 10 
y_upper <- max(plot_data$cum_mortality, na.rm = TRUE) * 1.1

# --- 主图 ---
p1 <- ggplot(plot_data, aes(x = time, y = cum_mortality, color = group)) +
  geom_step(linewidth = 1.2) +
  scale_color_manual(values = group_colors) +
  labs(x = "随访时间 (年)", y = "累积全因死亡率 (1/100,000)", color = "风险分组") +
  scale_x_continuous(breaks = target_times, limits = c(0, max_x)) +
  scale_y_continuous(limits = c(0, y_upper), expand = c(0,0)) +
  theme_minimal(base_family = "SimHei", base_size = 16) + 
  theme(
    legend.position = c(0.05, 0.95), 
    legend.justification = c(0, 1),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 20)
  ) +
  annotate("text", x = max_x, y = y_upper * 0.05, label = p_text, 
           hjust = 1, fontface = "bold", family = "SimHei", size = 5)

# 7. 绝对正确的风险表计算
km_summary <- summary(fit_km, times = target_times, extend = TRUE)
risk_table_data <- data.frame(
  time = km_summary$time,
  risk_group = gsub("risk_group=", "", as.character(km_summary$strata)),
  n = km_summary$n.risk
)

# 确保 risk_group 为因子（保持顺序）
risk_table_data$risk_group <- factor(risk_table_data$risk_group,
                                     levels = c("低危未筛查组", "高危筛查组", "高危未筛查组"))

# 绘制风险表
p2 <- ggplot(risk_table_data, aes(x = time, y = risk_group, label = n)) +
  geom_text(size = 5, family = "SimHei") +
  theme_void(base_family = "SimHei", base_size = 16) +
  scale_x_continuous(breaks = target_times, limits = c(0, max_x)) +
  labs(title = "人数统计") + 
  theme(
    plot.title = element_text(size = 16, face = "bold", family = "SimHei", hjust = 0),
    axis.text.y = element_text(hjust = 1, family = "SimHei"),
    plot.margin = margin(t = 0, r = 20, b = 20, l = 10) 
  )

# 8. 组合图片与保存 (⭐⭐ 核心修改：保存文件名变更 ⭐⭐)
final_plot <- p1 / p2 + plot_layout(heights = c(4, 1))

output_dir <- "D:/HuaweiMoveData/Users/Ai/Desktop/非吸烟小论文/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# 更新了文件名，标注 10年截断
output_path <- paste0(output_dir, "全因死亡率_KM版_10年截断.png")

ggsave(output_path, final_plot, width = 10, height = 8, dpi = 300)

cat(paste0("\n✅ 分析完成，图片已成功保存至：", output_path, "\n"))