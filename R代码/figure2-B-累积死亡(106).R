# ===================== 0. 环境彻底净化 =====================
rm(list = ls())
gc()

# ==============================================================================
# 1. 加载必需包
# ==============================================================================
library(cmprsk)   
library(dplyr)    
library(ggplot2)  
library(gridExtra)
library(tidyr)
library(lubridate)
library(data.table)
library(showtext) 

# ==============================================================================
# 2. 字体设置 (SimHei)
# ==============================================================================
showtext_auto()
font_add("SimHei", "simhei.ttf") # 系统中有simhei.ttf

# ==============================================================================
# 3. 数据读取与配置 (不再使用日期)
# ==============================================================================
input_file_path <- "D:/匹配结果_最终修正版.csv"

cat("正在读取数据...\n")
raw_data <- fread(input_file_path, header = TRUE, quote = "", fill = TRUE, data.table = FALSE, na.strings = c("", "NA", "NULL"))

# ==============================================================================
# 4. 数据清洗 (直接使用 death_time)
# ==============================================================================
data_clean <- raw_data %>%
  select(
    death_time = death_time,         # 直接读取随访天数
    lc_death   = LC_death,           # 肺癌死亡
    all_death  = all_cause_death,    # 全因死亡
    group      = as_lc_risk
  ) %>%
  mutate(
    # 将天数除以 365.25 换算成“年”
    time_years = as.numeric(as.character(death_time)) / 365.25,
    lc_death   = as.numeric(as.character(lc_death)),
    all_death  = as.numeric(as.character(all_death)),
    group      = as.numeric(as.character(group))
  ) %>%
  # 剔除随访时间缺失或无效的数据
  filter(!is.na(time_years) & time_years > 0)

# ==============================================================================
# 5. 定义竞争风险状态 (⭐⭐ 核心逻辑修改：增加 10 年截断 ⭐⭐)
# ==============================================================================
data_final <- data_clean %>%
  mutate(
    # 1. 原始状态判定
    status = case_when(
      lc_death == 1 ~ 1,                 # 状态1：肺癌死亡 (目标)
      all_death == 1 & lc_death == 0 ~ 2,# 状态2：死于其他原因 (竞争)
      TRUE ~ 0                           # 状态0：存活/删失
    )
  ) %>%
  mutate(
    # 2. 十年截断逻辑：时间大于10年的，状态转为删失，时间截断为10
    status = ifelse(time_years > 10, 0, status),
    time_years = ifelse(time_years > 10, 10, time_years)
  ) %>%
  filter(!is.na(group))

# 定义分组标签
data_final$risk_group <- factor(
  data_final$group,
  levels = c(0, 1, 2),
  labels = c("低危未筛查组", "高危筛查组", "高危未筛查组") 
)

cat("--- 事件分布 (10年截断后) ---\n")
print(table(data_final$status)) # 检查一下 1 和 2 的数量

# ==============================================================================
# 6. 计算 Fine-Gray 模型
# ==============================================================================
fit <- cuminc(ftime = data_final$time_years, fstatus = data_final$status, group = data_final$risk_group, cencode = 0)

# 提取数据
target_names <- names(fit)[sapply(names(fit), function(x) grepl(" 1$", x))]
cif_data <- data.frame()
for (nm in target_names) {
  group_name <- gsub(" 1$", "", nm)
  temp_df <- data.frame(time = fit[[nm]]$time, cif = fit[[nm]]$est * 100000, group = group_name)
  cif_data <- rbind(cif_data, temp_df)
}
cif_data$group <- factor(cif_data$group, levels = c("低危未筛查组", "高危筛查组", "高危未筛查组"))

# ==============================================================================
# 7. 绘图 
# ==============================================================================
p_value <- fit$Tests[1, "pv"]
p_text <- ifelse(p_value < 0.001, "Gray's 检验 P < 0.001", paste0("Gray's 检验 P = ", round(p_value, 3)))

group_colors <- c("低危未筛查组" = "#95e1d3", "高危筛查组" = "#f38181", "高危未筛查组" = "#fce38a")

# (⭐⭐ 修改点：强制 X 轴最大值为 10 ⭐⭐)
max_x <- 10
y_upper <- max(cif_data$cif, na.rm=TRUE) * 1.1

# --- 主图 ---
p1 <- ggplot(cif_data, aes(x = time, y = cif, color = group)) +
  geom_step(linewidth = 1.2) +
  scale_color_manual(values = group_colors) +
  labs(x = "随访时间 (年)", y = "累积死亡率 (1/100,000)", color = "风险分组") +
  scale_x_continuous(breaks = c(0, 2.5, 5, 7.5, 10), limits = c(0, max_x)) +
  scale_y_continuous(limits = c(0, y_upper), expand = c(0,0)) +
  theme_minimal(base_family = "SimHei") + 
  theme(
    legend.position = c(0.05, 0.95), 
    legend.justification = c(0, 1),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  # P值右下角
  annotate("text", 
           x = max_x,          
           y = y_upper * 0.05, 
           label = p_text, 
           hjust = 1,          
           vjust = 0,          
           fontface = "bold", 
           family = "SimHei")

# --- 风险表 ---
target_times <- c(0, 2.5, 5, 7.5, 10)
risk_table_data <- data.frame()

# 遍历每个目标时间点
for(t in target_times) {
  # 计算每个分组在时间t的风险人数：
  tmp <- data_final %>% 
    group_by(risk_group) %>% 
    summarise(
      n = sum(time_years >= t & !(status %in% c(1,2) & time_years < t)), 
      .groups = 'drop'
    ) %>% 
    mutate(time = t)
  risk_table_data <- rbind(risk_table_data, tmp)
}

p2 <- ggplot(risk_table_data, aes(x = time, y = risk_group, label = n)) +
  geom_text(size = 3.8, family = "SimHei") + 
  theme_void(base_family = "SimHei") + 
  scale_x_continuous(breaks = target_times, limits = c(0, max_x)) +
  labs(title = "  人数统计") + 
  theme(
    plot.title = element_text(size = 10, face = "bold", family = "SimHei"),
    axis.text.y = element_text(size = 10, hjust = 1, margin = margin(r = 10), family = "SimHei"),
    plot.margin = margin(t = 0, r = 10, b = 10, l = 10)
  )

# --- 组合输出 ---
final_plot <- grid.arrange(p1, p2, nrow = 2, heights = c(4, 1))

output_dir <- "D:/HuaweiMoveData/Users/Ai/Desktop/非吸烟小论文/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# 为避免覆盖原图，加上“十年截断”标识
output_path <- paste0(output_dir, "累积肺癌死亡率_竞争风险_十年截断_中文版.png")

# 保存图片
ggsave(output_path, final_plot, width = 8, height = 6, dpi = 300)
cat("\n✅ 分析完成！图片已成功保存至：", output_path, "\n")