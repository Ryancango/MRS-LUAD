# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/TOX_final.R
# Role: Analyze TOX and immune-checkpoint marker expression in CD8 T-cell subsets.
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
# 0. 环境与设置
# ==============================================================================
rm(list = ls())
setwd("~/lung_cancer/scRNA/131907")

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(cowplot)

print(paste0("[", Sys.time(), "] 正在加载数据..."))
seurat_obj <- readRDS("seurat_tlung_riskgroup.rds")

# ==============================================================================
# 1. 数据准备 (标准流程)
# ==============================================================================
# --- 恢复样本分组 ---
meta_epi <- subset(seurat_obj, subset = Cell_type == "Epithelial cells")@meta.data
sample_stats <- meta_epi %>%
  group_by(Sample) %>%
  summarise(n_high = sum(risk_group == "High Risk"), total = n(), ratio = n_high/total)

threshold <- median(sample_stats$ratio)
sample_stats$Sample_Group <- ifelse(sample_stats$ratio > threshold, "High Risk", "Low Risk")
risk_map <- setNames(sample_stats$Sample_Group, sample_stats$Sample)
seurat_obj@meta.data$Sample_Group <- unname(risk_map[as.character(seurat_obj@meta.data$Sample)])
seurat_obj <- subset(seurat_obj, subset = !is.na(Sample_Group))

# --- 提取 CD8 T 细胞 ---
cd8_cells <- subset(seurat_obj, subset = Cell_subtype %in% c("Exhausted CD8+ T", "CD8 low T", "Cytotoxic CD8+ T", "Naive CD8+ T"))
cd8_cells$Sample_Group <- factor(cd8_cells$Sample_Group, levels = c("Low Risk", "High Risk"))

# ==============================================================================
# 2. 左图：DotPlot (保留全貌)
# ==============================================================================
# 既然 HAVCR2 在样本级不显著，我们在 DotPlot 里依然保留它，
# 因为 DotPlot 展示的是"单细胞级"的趋势，那里它还是红的，能呼应 CellChat
genes_terminal <- c("TOX", "HAVCR2", "ENTPD1", "CTLA4") 
genes_progenitor <- c("TCF7", "LAG3", "TIGIT", "GZMK") 
genes_to_plot <- c(genes_terminal, genes_progenitor)

p_dot <- DotPlot(cd8_cells, 
                 features = genes_to_plot, 
                 group.by = "Sample_Group",
                 dot.scale = 8,
                 scale = TRUE) +  
  scale_color_gradientn(colors = c("#F0F0F0", "#FFD9D9", "#E64B35", "#8B0000"), name = "Avg Exp") +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey60", size = 0.5) +
  annotate("text", x = 2.5, y = 2.8, label = "Terminal", fontface = "bold", size = 4) +
  annotate("text", x = 6.5, y = 2.8, label = "Progenitor", fontface = "bold", size = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    axis.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    panel.grid.major = element_line(color = "grey95"), 
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(title = "Differential Exhaustion Phenotypes") +
  coord_cartesian(ylim = c(1, 3), clip = "off") 

# ==============================================================================
# 3. 右图：Boxplot (只画 TOX)
# ==============================================================================
# 准备数据
tox_exp <- AverageExpression(cd8_cells, features = c("TOX"), group.by = "Sample")$RNA
tox_df <- as.data.frame(t(tox_exp))
tox_df$Sample_Key <- rownames(tox_df)
tox_df$Sample <- gsub("-", "_", tox_df$Sample_Key)
tox_df$Sample <- gsub("\\.", "_", tox_df$Sample)

plot_data <- left_join(tox_df, sample_stats[, c("Sample", "Sample_Group")], by = "Sample")
plot_data <- plot_data %>% filter(!is.na(Sample_Group))
plot_data$Sample_Group <- factor(plot_data$Sample_Group, levels = c("Low Risk", "High Risk"))

# 绘图 (TOX)
p_tox <- ggplot(plot_data, aes(x = Sample_Group, y = TOX, fill = Sample_Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, size = 0.5, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.15, size = 3, shape = 21, fill = "white", color = "black", stroke = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x.npc = 0.5, vjust = 0.5, size = 6) +
  scale_fill_manual(values = c("High Risk" = "#E64B35", "Low Risk" = "#4DBBD5")) +
  theme_classic() +
  labs(title = "TOX Expression", y = "Avg Expression", x = "") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.text.x = element_text(size = 11, color = "black", face="bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.line = element_line(size = 0.5),
    legend.position = "none"
  )

# ==============================================================================
# 4. 拼图 (左右结构)
# ==============================================================================
pdf("Fig7_TOX_Final_Clean.pdf", width = 10, height = 6)

# 布局：左图占 1.8，右图占 1
layout <- p_dot + p_tox + plot_layout(widths = c(1.8, 1))

print(layout)
dev.off()

print("✅ 精简版耗竭图已生成：Fig7_TOX_Final_Clean.pdf")




# ==============================================================================
# 1. 数据准备
# ==============================================================================
# --- 恢复样本分组 ---
meta_epi <- subset(seurat_obj, subset = Cell_type == "Epithelial cells")@meta.data
sample_stats <- meta_epi %>%
  group_by(Sample) %>%
  summarise(n_high = sum(risk_group == "High Risk"), total = n(), ratio = n_high/total)

threshold <- median(sample_stats$ratio)
sample_stats$Sample_Group <- ifelse(sample_stats$ratio > threshold, "High Risk", "Low Risk")
risk_map <- setNames(sample_stats$Sample_Group, sample_stats$Sample)
seurat_obj@meta.data$Sample_Group <- unname(risk_map[as.character(seurat_obj@meta.data$Sample)])
seurat_obj <- subset(seurat_obj, subset = !is.na(Sample_Group))

# --- 提取 CD8 T 细胞 ---
cd8_cells <- subset(seurat_obj, subset = Cell_subtype %in% c("Exhausted CD8+ T", "CD8 low T", "Cytotoxic CD8+ T", "Naive CD8+ T"))
cd8_cells$Sample_Group <- factor(cd8_cells$Sample_Group, levels = c("Low Risk", "High Risk"))

# ==============================================================================
# 2. 左图：DotPlot (TOX 放最后 + 顶部箭头)
# ==============================================================================
# --- 1. 基因排序：初始 -> 中间 -> 终末 -> TOX(压轴) ---
genes_progenitor <- c("TCF7", "IL7R")
genes_intermediate <- c("GZMK", "LAG3", "TIGIT")
# 将 TOX 放在 Terminal 组的最后面
genes_terminal <- c("PDCD1", "HAVCR2", "ENTPD1", "CTLA4", "TOX")

genes_to_plot <- c(genes_progenitor, genes_intermediate, genes_terminal)

# --- 2. 绘图 ---
p_dot <- DotPlot(cd8_cells, 
                 features = genes_to_plot, 
                 group.by = "Sample_Group",
                 dot.scale = 8,
                 scale = TRUE) +  
  scale_color_gradientn(colors = c("#F0F0F0", "#FFD9D9", "#E64B35", "#8B0000"), name = "Avg Exp") +
  
  # --- 顶部箭头与文字 ---
  # 箭头从第1个基因指向最后一个基因
  annotate("segment", x = 1, xend = length(genes_to_plot), y = 3.2, yend = 3.2,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"), 
           color = "black", size = 0.8) +
  # 文字标签
  annotate("text", x = 1, y = 3.5, label = "Progenitor", hjust = 0.2, fontface = "bold", size = 4) +
  annotate("text", x = length(genes_to_plot), y = 3.5, label = "Terminal Exhaustion", hjust = 0.8, fontface = "bold", size = 4) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 11, face = "italic"),
    axis.text.y = element_text(color = "black", size = 12, face = "bold"),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    # 这里先不强制设置 legend 位置，交给 patchwork 统一管理
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10) 
  ) +
  labs(title = "Exhaustion Trajectory") +
  coord_cartesian(ylim = c(1, 3), clip = "off")

# ==============================================================================
# 3. 右图：2x2 Boxplot (保持不变)
# ==============================================================================
target_genes <- c("TOX", "PDCD1", "HAVCR2", "ENTPD1") 

avg_exp <- AverageExpression(cd8_cells, features = target_genes, group.by = "Sample")$RNA
plot_df <- as.data.frame(t(avg_exp))
plot_df$Sample_Key <- rownames(plot_df)
plot_df$Sample <- gsub("-", "_", plot_df$Sample_Key)
plot_df$Sample <- gsub("\\.", "_", plot_df$Sample)

box_data <- left_join(plot_df, sample_stats[, c("Sample", "Sample_Group")], by = "Sample")
box_data <- box_data %>% filter(!is.na(Sample_Group))
box_data$Sample_Group <- factor(box_data$Sample_Group, levels = c("Low Risk", "High Risk"))

create_boxplot <- function(data, gene, title_label) {
  if(!gene %in% colnames(data)) return(ggplot() + theme_void())
  ggplot(data, aes_string(x = "Sample_Group", y = gene, fill = "Sample_Group")) +
    geom_boxplot(width = 0.5, outlier.shape = NA, size = 0.4, alpha = 0.8, color = "black") +
    geom_jitter(width = 0.15, size = 2, shape = 21, fill = "white", color = "black", stroke = 0.3) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", label.x.npc = 0.5, vjust = 0.5, size = 5) +
    scale_fill_manual(values = c("High Risk" = "#E64B35", "Low Risk" = "#4DBBD5")) +
    theme_classic() +
    labs(title = title_label, y = "Avg Exp", x = "") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 12),
      axis.text.x = element_blank(), 
      axis.text.y = element_text(size = 9, color = "black"),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 9),
      legend.position = "none", # 确保箱线图没有图例干扰
      plot.margin = margin(5, 5, 5, 5)
    )
}

p1 <- create_boxplot(box_data, "TOX", "TOX") + theme(axis.text.x = element_blank())
p2 <- create_boxplot(box_data, "PDCD1", "PDCD1 (PD-1)") + theme(axis.text.x = element_blank())
p3 <- create_boxplot(box_data, "HAVCR2", "HAVCR2 (TIM-3)") + theme(axis.text.x = element_text(size=10, color="black", face="bold"))
p4 <- create_boxplot(box_data, "ENTPD1", "ENTPD1 (CD39)") + theme(axis.text.x = element_text(size=10, color="black", face="bold"))

p_right_grid <- (p1 + p2) / (p3 + p4)

# ==============================================================================
# 4. 最终拼图 (修复图例显示)
# ==============================================================================
# 增加高度以容纳底部图例
pdf("Fig7_Final_TOX_End.pdf", width = 11, height = 7)

# 使用 guides = 'collect' 强制将图例放在最底部
# theme(legend.position = "bottom") 确保它在下方
layout <- p_dot + p_right_grid + 
  plot_layout(widths = c(1.4, 1), guides = 'collect') & 
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.margin = margin(t = 10))

print(layout)
dev.off()

print("✅ 完美版图表已生成：Fig7_Final_TOX_End.pdf")
