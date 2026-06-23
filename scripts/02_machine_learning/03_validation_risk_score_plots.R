# Curated release script for: MRS-LUAD
# Source script: intersection_mechine/mechine/outside_data/验证集风险评分.R
# Role: Generate validation risk-score distributions and survival status plots.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}

# ==============================================================================
# 1. 紧急数据恢复 (根据你之前的代码逻辑复原)
# ==============================================================================
library(survival)
library(dplyr)

# 你的 10 个核心基因 (用于 filter_data)
final_genes <- c("RRAS", "GPM6B", "FMO3", "CHRDL1", "ADRB2", 
                 "LPL", "ANGPT1", "SPOCK2", "P2RY14", "ADH1B")

# 定义过滤函数
filter_data <- function(df, genes) {
  # 确保包含生存列
  survival_cols <- c("OS.time", "OS")
  # 确保包含基因列
  gene_cols <- intersect(genes, colnames(df))
  
  # 移除 OS.time <= 0
  df <- df[df$OS.time > 0, ]
  
  # 返回子集
  return(df[, c(survival_cols, gene_cols)])
}

# 过滤数据 (假设原始大对象还在，如果不在请先 read.csv 读取)
# 如果这三个对象也没了，请先重新读取原始csv/txt文件
if(exists("tcgaluad") & exists("GSE31210") & exists("GSE72094")) {
  
  tcga_filtered <- filter_data(tcgaluad, final_genes)
  gse31210_filtered <- filter_data(GSE31210, final_genes)
  gse72094_filtered <- filter_data(GSE72094, final_genes)
  
  # --- 关键：使用相同的随机种子复原拆分 ---
  set.seed(5201314)
  
  # 1. TCGA Split
  train_idx_tcga <- sample(nrow(tcga_filtered), size = round(0.5 * nrow(tcga_filtered)))
  tcga_train <- tcga_filtered[train_idx_tcga, ]
  tcga_vali <- tcga_filtered[-train_idx_tcga, ]
  
  # 2. GSE31210 Split
  train_idx_gse31210 <- sample(nrow(gse31210_filtered), size = round(0.5 * nrow(gse31210_filtered)))
  gse31210_train <- gse31210_filtered[train_idx_gse31210, ]
  gse31210_vali <- gse31210_filtered[-train_idx_gse31210, ]
  
  # 3. GSE72094 Split
  train_idx_gse72094 <- sample(nrow(gse72094_filtered), size = round(0.5 * nrow(gse72094_filtered)))
  gse72094_train <- gse72094_filtered[train_idx_gse72094, ]
  gse72094_vali <- gse72094_filtered[-train_idx_gse72094, ]
  
  # 4. 重建训练集
  train_data <- rbind(tcga_train, gse31210_train, gse72094_train)
  
  print(paste("✅ 数据恢复成功！Train Data 样本数:", nrow(train_data)))
  
} else {
  stop("❌ 错误：原始数据集 (tcgaluad, GSE31210, GSE72094) 丢失。请先重新读取这三个原始文件！")
}

# ==============================================================================
# 2. 重建 Cox 模型
# ==============================================================================
cox_formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(final_genes, collapse = "+")))
final_cox_model <- coxph(cox_formula, data = train_data)
print("✅ 模型重建成功！")



# ==============================================================================
# 0. 环境与数据准备
# ==============================================================================
library(ggplot2)
library(dplyr)
library(scales)

# 假设 plot_data 已经生成 (如果不在内存中，请运行下面的快速生成代码)
# --- 快速生成 plot_data (以防万一) ---
if(!exists("merged_vali")) {
  # (这里省略数据恢复步骤，假设你已经运行了上一条回复的代码)
  stop("请先运行上一条代码恢复 merged_vali 数据！")
}
# 重新计算 plot_data
plot_data <- merged_vali %>% arrange(RiskScore)
plot_data$Rank <- 1:nrow(plot_data)
cutoff <- median(plot_data$RiskScore)
plot_data$RiskGroup <- ifelse(plot_data$RiskScore > cutoff, "High Risk", "Low Risk")
plot_data$Status <- ifelse(plot_data$OS == 1, "Dead", "Alive")
# -------------------------------------

# ==============================================================================
# 1. 计算双轴转换系数 (Scaling Factor)
# ==============================================================================
# 目的：把 Survival Time (0~3000) 压缩到 Risk Score (-5~5) 的范围内
# 公式：y_scaled = (y - min_y) / (max_y - min_y) * (max_x - min_x) + min_x
# 为了简单且美观，我们使用线性比例缩放

# 获取两个变量的范围
range_risk <- range(plot_data$RiskScore)
range_time <- range(plot_data$OS.time)

# 计算缩放系数 (scale) 和 偏移量 (shift)
# 我们希望 Time 的 0 点对齐 Risk 的最小值，Time 的最大值对齐 Risk 的最大值
scale_factor <- (range_risk[2] - range_risk[1]) / (range_time[2] - range_time[1])
shift_value <- range_risk[1] - range_time[1] * scale_factor

# 定义转换函数 (用于画图)
trans_time_to_risk <- function(x) { x * scale_factor + shift_value }

# 定义逆转换函数 (用于坐标轴标签)
trans_risk_to_time <- function(x) { (x - shift_value) / scale_factor }

# ==============================================================================
# 2. 绘制双轴图
# ==============================================================================
# 配色定义
col_risk <- c("Low Risk" = "#4DBBD5", "High Risk" = "#E64B35") # 风险分颜色
col_status <- c("Alive" = "grey70", "Dead" = "black")          # 生存状态颜色 (用黑灰区分，避免颜色过多杂乱)

p_dual <- ggplot(plot_data, aes(x = Rank)) +
  
  # --- 图层 1: 生存时间 (右轴数据) ---
  # 注意：必须先对 y 进行转换，使其落入左轴范围
  # 放在底层，用竖线或点表示
  geom_segment(aes(xend = Rank, y = min(RiskScore), yend = trans_time_to_risk(OS.time), color = Status), 
               alpha = 0.3, size = 0.2) + # 加一条淡淡的垂线，增加高级感
  geom_point(aes(y = trans_time_to_risk(OS.time), shape = Status, color = Status), 
             size = 1.5, alpha = 0.6) +
  
  # --- 图层 2: 风险评分 (左轴数据) ---
  # 放在上层，因为是主要排序依据
  geom_point(aes(y = RiskScore, fill = RiskGroup), 
             shape = 21, color = "white", size = 2.5, stroke = 0.2) +
  
  # --- 辅助线: Cutoff ---
  geom_vline(xintercept = sum(plot_data$RiskGroup == "Low Risk") + 0.5, 
             linetype = "dashed", color = "grey40", size = 0.8) +
  
  # --- 坐标轴设置 (核心) ---
  scale_y_continuous(
    # 左轴名称
    name = "Risk Score (S-Curve)",
    # 右轴配置：定义变换关系
    sec.axis = sec_axis(trans = ~ trans_risk_to_time(.), 
                        name = "Overall Survival Time (Days)",
                        breaks = pretty_breaks(n = 5))
  ) +
  
  # --- 颜色与形状映射 ---
  scale_fill_manual(values = col_risk, name = "Risk Group") + # 对应 RiskScore (实心点填充)
  scale_color_manual(values = col_status, name = "Survival Status") + # 对应 Time (边框/垂线)
  scale_shape_manual(values = c("Alive" = 1, "Dead" = 16), name = "Survival Status") + # 空心圆活，实心圆死
  
  # --- 主题美化 ---
  theme_classic() +
  labs(title = "Integrated Risk & Survival Landscape", x = "Patients (Ranked by Risk Score)") +
  
  theme(
    # 标题
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    
    # 左轴 (Risk Score) 样式 -> 对应红蓝色
    axis.title.y.left = element_text(face = "bold", size = 12, color = "#E64B35"),
    axis.text.y.left = element_text(color = "#E64B35", size = 10),
    axis.line.y.left = element_line(color = "#E64B35"),
    
    # 右轴 (Survival Time) 样式 -> 对应黑灰色
    axis.title.y.right = element_text(face = "bold", size = 12, color = "black", angle = 90),
    axis.text.y.right = element_text(color = "black", size = 10),
    axis.line.y.right = element_line(color = "black"),
    
    # X轴
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # 图例
    legend.position = "top",
    legend.box = "horizontal",
    legend.background = element_rect(fill = "white", color = "black", size = 0.2)
  )

# ==============================================================================
# 3. 拼图 (Dual Plot + Heatmap)
# ==============================================================================
# 如果你还想保留热图在下面，可以拼在一起

# (重新生成热图 p3，代码同上一次)
exp_mat <- scale(plot_data[, final_genes])
exp_mat[exp_mat > 2] <- 2; exp_mat[exp_mat < -2] <- -2
exp_df <- as.data.frame(exp_mat); exp_df$Rank <- plot_data$Rank
exp_long <- pivot_longer(exp_df, cols = all_of(final_genes), names_to = "Gene", values_to = "Expression")
exp_long$Gene <- factor(exp_long$Gene, levels = final_genes)

p3 <- ggplot(exp_long, aes(x = Rank, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_bw() +
  labs(x = "Patients (Ranked)", y = "") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 9, face = "bold.italic", color = "black"),
        panel.grid = element_blank(), legend.position = "right")

# 拼接：上图是双轴图，下图是热图
final_dual_plot <- p_dual / p3 + plot_layout(heights = c(2, 1.2))

# 保存
pdf("Fig_Merged_Validation_Dual_Axis.pdf", width = 10, height = 10)
print(final_dual_plot)
dev.off()

print("✅ 双Y轴合并图已生成：Fig_Merged_Validation_Dual_Axis.pdf")