# ===================== 0. 彻底净化环境与加载包 =====================
rm(list = ls())
gc()
cat("\n✨ 环境已彻底净化，开始运行 7 变量 IPW 校正代码...\n")

# 加载包
library(tableone)
library(survey)
library(dplyr)
library(tidyr)
library(data.table)

# ===================== 1. 读取数据 =====================
file_path <- "D:/匹配结果_最终修正版.csv"

cat("\n🚀 正在尝试使用 fread 高效读取...\n")

data <- tryCatch({
  fread(file_path, 
        encoding = "GBK", 
        quote = "",      # 禁用引号解析，防断行
        fill = TRUE, 
        check.names = FALSE,
        data.table = FALSE) 
}, error = function(e) {
  cat("GBK读取失败，尝试UTF-8...\n")
  fread(file_path, 
        encoding = "UTF-8", 
        quote = "", 
        fill = TRUE, 
        check.names = FALSE,
        data.table = FALSE)
})

cat("✅ 读取完成，当前总行数：", nrow(data), "\n")

# ===================== 2. 定义目标变量与预处理 =====================
# 【核心修改】：只保留你需要的 7 个变量
baseline_vars <- c(
  "gender",                          # 性别
  "age_group",                       # 年龄分组
  "all_second_smoke",                # 二手烟暴露
  "as_occupational_hazard_exposure", # 职业暴露
  "as_chronic_respiratory_disease",  # 慢性呼吸疾病史
  "yijijiazushi",                    # 一级家族史
  "BMI_GROUP"                        # BMI分组
)

# 筛选 as_lc_risk 为 1 和 2 的人群
data$as_lc_risk <- as.numeric(as.character(data$as_lc_risk))
data_filtered <- data %>% filter(as_lc_risk %in% c(1, 2))

# 重编码：1->1 (筛查), 2->0 (未筛查)
data_filtered$as_lc_risk <- ifelse(data_filtered$as_lc_risk == 1, 1, 0)

# BMI计算与分组 (如果原本没有 BMI_GROUP 列，在这里生成)
data_filtered <- data_filtered %>%
  mutate(
    BMI = as.numeric(as.character(BMI)), 
    BMI_GROUP = case_when(
      BMI < 18.5 ~ "<18.5",                  
      BMI >= 18.5 & BMI < 24.0 ~ "18.5-24.0",
      BMI >= 24.0 & BMI < 28.0 ~ "24.0-28.0",
      BMI >= 28.0 ~ "≥28.0",                 
      TRUE ~ NA_character_                   
    ) %>% factor(levels = c("<18.5", "18.5-24.0", "24.0-28.0", "≥28.0"))
  )

# ===================== 3. 缺失值剔除与格式转换 =====================
ps_cols <- c("as_lc_risk", baseline_vars)

# 确保 7 个基线变量全是 Factor 类型 (TableOne 计算分类变量 SMD 的关键)
data_filtered[baseline_vars] <- lapply(data_filtered[baseline_vars], as.factor)

# 剔除这 7 个变量的缺失值
data_complete <- data_filtered %>% drop_na(all_of(ps_cols))

cat("\n🎯 剔除缺失值后参与 IPW 的最终人数：", nrow(data_complete), "\n")
cat("  - 筛查组 (1):", sum(data_complete$as_lc_risk == 1), "\n")
cat("  - 未筛查组 (0):", sum(data_complete$as_lc_risk == 0), "\n")

# ===================== 4. IPW 计算 =====================
# 只针对这 7 个变量拟合倾向性评分模型
formula_ps <- as.formula(paste("as_lc_risk ~", paste(baseline_vars, collapse = " + ")))
ps_model <- glm(formula_ps, data = data_complete, family = binomial())

data_complete$ps <- predict(ps_model, type = "response")

# 计算 ATE 倒数权重
data_complete$ipw <- ifelse(data_complete$as_lc_risk == 1, 
                            1 / data_complete$ps, 
                            1 / (1 - data_complete$ps))

# 权重截断 (1% - 99%)，防止极端权重放大方差
q_low <- quantile(data_complete$ipw, 0.01)
q_high <- quantile(data_complete$ipw, 0.99)
data_complete$ipw_truncated <- pmin(pmax(data_complete$ipw, q_low), q_high)

# ===================== 5. TableOne 计算与结果提取 =====================
cat("\n正在计算加权前后的基线表与 SMD...\n")

# 5.1 加权前 (Unadjusted)
tab_unadj <- CreateTableOne(vars = baseline_vars, 
                            strata = "as_lc_risk", 
                            data = data_complete, 
                            test = FALSE) 

# 5.2 加权后 (Adjusted IPW)
svy_design <- svydesign(ids = ~1, data = data_complete, weights = ~ipw_truncated)
tab_adj <- svyCreateTableOne(vars = baseline_vars, 
                             strata = "as_lc_risk", 
                             data = svy_design, 
                             test = FALSE)

# 提取数据矩阵
mat_unadj <- print(tab_unadj, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE)
mat_adj   <- print(tab_adj, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE)

# 提取各部分数值
smd_unadj <- as.numeric(mat_unadj[, "SMD"])  
smd_adj   <- as.numeric(mat_adj[, "SMD"])    
unadj_0   <- mat_unadj[, "0"]                
unadj_1   <- mat_unadj[, "1"]                
adj_0     <- mat_adj[, "0"]                  
adj_1     <- mat_adj[, "1"]                  

# 构建最终输出的 Data Frame
final_df <- data.frame(
  Feature = rownames(mat_unadj),
  Level   = mat_unadj[, "level"],
  Unscreened_Unadj = unadj_0,
  Screened_Unadj   = unadj_1,
  Unscreened_Adj   = adj_0,
  Screened_Adj     = adj_1,
  SMD_Unadjusted   = sprintf("%.3f", smd_unadj),
  SMD_Adjusted     = sprintf("%.3f", smd_adj)
)

# 清理格式
rownames(final_df) <- NULL
final_df$Level[is.na(final_df$Level)] <- ""

# 中文化列名（自带最终样本量）
colnames(final_df) <- c(
  "变量", 
  "水平", 
  paste0("加权前未筛查组 (N=", sum(data_complete$as_lc_risk==0), ")"),
  paste0("加权前筛查组 (N=", sum(data_complete$as_lc_risk==1), ")"),
  "加权后未筛查组",
  "加权后筛查组",
  "加权前SMD", 
  "加权后SMD"
)

# ===================== 6. 导出最终结果 =====================
# 将输出文件名改为明确的 7 变量版本
out_path <- "D:/7变量校正_IPW基线表.csv"
write.csv(final_df, out_path, row.names = FALSE, fileEncoding = "GBK") 

cat("\n🎉 分析圆满完成！\n")
cat("📁 结果已成功导出至：", out_path, "\n")