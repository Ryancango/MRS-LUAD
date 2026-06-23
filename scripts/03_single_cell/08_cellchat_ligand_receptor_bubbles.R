# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/CellChat_Analysis/cellchatplot_bubble.R
# Role: Draw CellChat ligand-receptor bubble plots.
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
# 0. 环境与数据加载
# ==============================================================================
rm(list = ls())
library(CellChat)
library(ggplot2)
library(dplyr)
library(stringr) # 用于提取配体名字进行排序

# 设置路径
work_dir <- "~/lung_cancer/scRNA/131907/CellChat_Analysis"
setwd(work_dir)

print(paste0("[", Sys.time(), "] 正在加载数据..."))
cellchat <- readRDS("cellchat_object_computed.rds")

# ==============================================================================
# 1. 提取通讯数据
# ==============================================================================
# 定义分组索引
group.new <- levels(cellchat@idents)
src_low <- which(group.new == "Epi_Low")
src_high <- which(group.new == "Epi_High")
# 定义接收细胞
target_names <- c("Endothelial cells", "Myeloid cells", "T lymphocytes", "Fibroblasts")
target_idx <- which(group.new %in% target_names)

print("正在提取显著的配受体对...")

# --- A. 提取 Low Risk 数据 ---
df_low <- subsetCommunication(cellchat, sources.use = src_low, targets.use = target_idx)
if(nrow(df_low) > 0) df_low$Risk <- "Low Risk"

# --- B. 提取 High Risk 数据 ---
df_high <- subsetCommunication(cellchat, sources.use = src_high, targets.use = target_idx)
if(nrow(df_high) > 0) df_high$Risk <- "High Risk"

# --- C. 合并数据 ---
plot_data <- rbind(df_low, df_high)

# ==============================================================================
# 2. 数据清洗与智能排序 (关键步骤)
# ==============================================================================
# 这里的目的是让同一个家族的配体（如VEGFA, VEGFB）排在一起，像你图里那样

# 1. 提取配体主名 (例如从 "VEGFA - VEGFR1" 中提取 "VEGFA")
plot_data <- plot_data %>%
  mutate(Ligand_Family = str_extract(interaction_name_2, "^[A-Za-z0-9]+"))

# 2. 筛选 Top 互作 (防止图太长)
# 逻辑：取概率最高的互作，或者你可以根据图里的内容手动指定 filter
# 这里我们保留所有提取到的数据，但如果太多可以取消下面两行的注释
# top_interactions <- plot_data %>% group_by(interaction_name_2) %>% summarize(m=max(prob)) %>% arrange(desc(m)) %>% head(50) %>% pull(interaction_name_2)
# plot_data <- plot_data %>% filter(interaction_name_2 %in% top_interactions)

# 3. 设置因子水平 (Factor Levels) 用于排序
# 我们先按 Ligand_Family 排序，再按 interaction_name_2 排序
# unique() 保持出现的顺序
sorted_interactions <- plot_data %>%
  arrange(Ligand_Family, interaction_name_2) %>%
  pull(interaction_name_2) %>%
  unique()

plot_data$interaction_name_2 <- factor(plot_data$interaction_name_2, levels = sorted_interactions)

# 4. 设置 Risk 分组顺序 (左 Low 右 High)
plot_data$Risk <- factor(plot_data$Risk, levels = c("Low Risk", "High Risk"))

# 5. 设置 Target 细胞顺序
plot_data$target <- factor(plot_data$target, levels = target_names)

# 6. 处理 P 值 (用于气泡大小)
# 避免 log(0) 为无穷大，加一个极小值
plot_data$p_val_log <- -log10(plot_data$pval + 1e-10)

# ==============================================================================
# 3. 绘图 (完全复刻样式)
# ==============================================================================
print("正在绘制高清气泡图...")

p <- ggplot(plot_data, aes(x = target, y = interaction_name_2)) +
  
  # 画气泡：
  # shape=21 是空心圆，可以分别设置填充色(fill)和边框色(color)
  # stroke=0.2 让边框细一点，显得精致
  geom_point(aes(size = p_val_log, fill = prob), shape = 21, color = "black", stroke = 0.2) +
  
  # 分面：左右对比
  facet_wrap(~Risk, ncol = 2) +
  
  # 配色方案：Spectral (蓝-黄-红)
  # 这种配色对比度高，非常适合展示强度
  scale_fill_gradientn(colors = c("#3288BD", "#99D594", "#E6F598", "#FEE08B", "#FC8D59", "#D53E4F"),
                       name = "Interaction\nStrength") +
  
  # 大小图例
  scale_size_continuous(name = "-log10(P-value)", range = c(2, 6)) +
  
  # 主题设置
  theme_bw() +
  theme(
    # 字体
    text = element_text(size = 10, color = "black"),
    
    # 坐标轴
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face="bold"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_blank(), # 去掉 X/Y 轴标题
    
    # 分面标题 (Low Risk / High Risk)
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    
    # 网格线 (极淡的灰色，辅助对齐)
    panel.grid.major = element_line(color = "grey95", size = 0.2),
    panel.grid.minor = element_blank(),
    
    # 图例
    legend.position = "right"
  ) +
  
  labs(title = "Differential Communication: Epi Outgoing Signaling")

# ==============================================================================
# 4. 保存
# ==============================================================================
# 自动计算高度，防止重叠
plot_height <- max(6, length(unique(plot_data$interaction_name_2)) * 0.28 + 1.5)

pdf("Fig6_Bubble_Reproduced_Spectral.pdf", width = 12, height = plot_height)
print(p)
dev.off()

print(paste0("✅ 图片已生成: Fig6_Bubble_Reproduced_Spectral.pdf (高度: ", round(plot_height, 1), " inch)"))