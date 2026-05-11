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
font_add("SimHei", "simhei.ttf") # 确保系统中有simhei.ttf，Windows默认有

# ==============================================================================
# 3. 数据读取与清洗
# ==============================================================================
input_file_path <- "D:/匹配结果_最终修正版.csv"

cat("正在读取数据...\n")
raw_data <- fread(input_file_path, header = TRUE, quote = "", fill = TRUE, data.table = FALSE, na.strings = c("", "NA", "NULL"))

data_clean <- raw_data %>%
  select(
    cancer_time = cancer_time,
    death_time  = death_time,
    lc_event    = LC_event,
    death       = all_cause_death,
    group       = as_lc_risk
  ) %>%
  mutate(
    # 将字符型转为数值型，并将天数除以 365.25 换算成“年”
    cancer_years = as.numeric(as.character(cancer_time)) / 365.25,
    death_years  = as.numeric(as.character(death_time)) / 365.25,
    lc_event     = as.numeric(as.character(lc_event)),
    death        = as.numeric(as.character(death)),
    group        = as.numeric(as.character(group))
  )

# ==============================================================================
# 4. 竞争风险模型状态定义与计算 (⭐⭐ 核心修改：10年截断逻辑 ⭐⭐)
# ==============================================================================
data_final <- data_clean %>%
  mutate(
    status = case_when(
      lc_event == 1 ~ 1,                # 状态1：发生肺癌发病（目标事件）
      death == 1 & lc_event == 0 ~ 2,   # 状态2：未得肺癌但全因死亡（竞争事件）
      TRUE ~ 0                          # 状态0：删失（存活且未发病）
    ),
    # 动态匹配随访时间：发生发病看cancer_years，发生死亡看death_years，删失默认取发病随访时间
    time_years = case_when(
      status == 1 ~ cancer_years,
      status == 2 ~ death_years,
      TRUE ~ coalesce(cancer_years, death_years) 
    )
  ) %>%
  # -------------------- 新增：十年截断处理 --------------------
# 如果实际随访时间大于 10 年，则在第 10 年时将其状态转为删失，时间固定为 10
mutate(
  status = ifelse(time_years > 10, 0, status),
  time_years = ifelse(time_years > 10, 10, time_years)
) %>%
  # ------------------------------------------------------------
# 剔除时间无效或分组缺失的人员
filter(!is.na(time_years), time_years > 0, !is.na(group))

# 因子化风险分组
data_final$risk_group <- factor(
  data_final$group,
  levels = c(0, 1, 2),
  labels = c("低危未筛查组", "高危筛查组", "高危未筛查组") 
)

cat("数据清洗完成，进入模型计算，当前样本量：", nrow(data_final), "\n")

# 计算累积发生率 (CIF)
fit <- cuminc(ftime = data_final$time_years, fstatus = data_final$status, group = data_final$risk_group, cencode = 0)

# 提取数据 (提取目标事件 status 1 的数据)
target_names <- names(fit)[sapply(names(fit), function(x) grepl(" 1$", x))]
cif_data <- data.frame()
for (nm in target_names) {
  group_name <- gsub(" 1$", "", nm)
  temp_df <- data.frame(time = fit[[nm]]$time, cif = fit[[nm]]$est * 100000, group = group_name)
  cif_data <- rbind(cif_data, temp_df)
}
cif_data$group <- factor(cif_data$group, levels = c("低危未筛查组", "高危筛查组", "高危未筛查组"))

# ==============================================================================
# 5. 绘图 
# ==============================================================================
p_value <- fit$Tests[1, "pv"]
p_text <- ifelse(p_value < 0.001, "Gray's 检验 P < 0.001", paste0("Gray's 检验 P = ", round(p_value, 3)))

group_colors <- c("低危未筛查组" = "#95e1d3", "高危筛查组" = "#f38181", "高危未筛查组" = "#fce38a")

# (⭐⭐ 修改点：强制最大X轴为 10 ⭐⭐)
max_x <- 10 
y_upper <- max(cif_data$cif, na.rm=TRUE) * 1.1

# --- 主图 ---
p1 <- ggplot(cif_data, aes(x = time, y = cif, color = group)) +
  geom_step(linewidth = 1.2) +
  scale_color_manual(values = group_colors) +
  labs(x = "随访时间 (年)", y = "累积发生率 (1/100,000)", color = "风险分组") +
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
for(t in target_times) {
  tmp <- data_final %>% 
    group_by(risk_group) %>% 
    summarise(n = sum(time_years >= t & !(status %in% c(1,2) & time_years < t)), 
              .groups = 'drop') %>% 
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

# 保存图片
output_dir <- "D:/HuaweiMoveData/Users/Ai/Desktop/非吸烟小论文/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_path <- paste0(output_dir, "累积肺癌发病率_竞争风险_十年截断_中文版.png")

ggsave(output_path, final_plot, width = 8, height = 6, dpi = 300)
cat("\n✅ 分析完成，图片已成功保存至：", output_path, "\n")