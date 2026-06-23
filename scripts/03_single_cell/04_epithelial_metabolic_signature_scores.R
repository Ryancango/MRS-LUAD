# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/代谢流hallmark箱线图.R
# Role: Compare epithelial-cell metabolic signatures between high- and low-MRS groups.
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
setwd("~/lung_cancer/scRNA/131907")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)

print(paste0("[", Sys.time(), "] 正在加载数据..."))
seurat_obj <- readRDS("seurat_tlung_riskgroup.rds")
seurat_epi <- subset(seurat_obj, subset = Cell_type == "Epithelial cells")

# ==============================================================================
# 1. 定义 4 个基因集并计算评分
# ==============================================================================
print("正在计算 4 大代谢/环境评分...")

# A. 缺氧 (Hypoxia) - 新增
genes_hypoxia <- c("VEGFA", "SLC2A1", "PGK1", "CA9", "LDHA", "ENO1", "ADM", "BNIP3", "NDRG1")

# B. 糖酵解 (Glycolysis)
genes_glycolysis <- c("LDHA", "ENO1", "HK2", "PKM", "PGK1", "GAPDH", "SLC2A1")

# C. 脂质合成 (Lipid Biosynthesis)
genes_lipid_syn  <- c("FASN", "SCD", "ACLY", "ACACA", "IDH1")

# D. 脂质摄取 (Lipid Uptake)
genes_lipid_uptake <- c("LPL", "CD36", "LIPE", "FABP4")

# 执行打分
seurat_epi <- AddModuleScore(seurat_epi, features = list(genes_hypoxia), name = "Score_Hypoxia")
seurat_epi <- AddModuleScore(seurat_epi, features = list(genes_glycolysis), name = "Score_Glyco")
seurat_epi <- AddModuleScore(seurat_epi, features = list(genes_lipid_syn), name = "Score_LipidSyn")
seurat_epi <- AddModuleScore(seurat_epi, features = list(genes_lipid_uptake), name = "Score_LipidUptake")

# 提取数据
plot_data <- seurat_epi@meta.data %>%
  select(risk_group, Score_Hypoxia1, Score_Glyco1, Score_LipidSyn1, Score_LipidUptake1)

colnames(plot_data) <- c("RiskGroup", "Hypoxia", "Glycolysis", "Lipid_Synthesis", "Lipid_Uptake")
plot_data$RiskGroup <- factor(plot_data$RiskGroup, levels = c("High Risk", "Low Risk"))

# ==============================================================================
# 2. 定义手动绘图函数 (保持之前的样式)
# ==============================================================================
my_cols <- c("High Risk" = "#E64B35", "Low Risk" = "#4DBBD5")

plot_box_strict_manual <- function(data, y_col, title_text, star_y, view_min, view_max) {
  
  ggplot(data, aes_string(x = "RiskGroup", y = y_col, fill = "RiskGroup")) +
    geom_boxplot(width = 0.5, outlier.shape = NA, size = 0.6, color = "black") +
    stat_compare_means(method = "wilcox.test", 
                       label = "p.signif", 
                       comparisons = list(c("Low Risk", "High Risk")),
                       label.y = star_y,
                       tip.length = 0.01,
                       size = 6) +
    scale_fill_manual(values = my_cols) +
    theme_classic() +
    coord_cartesian(ylim = c(view_min, view_max), clip = "off") +
    labs(title = title_text, y = "Signature Score", x = "") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "none",
      plot.margin = margin(t = 30, r = 10, b = 5, l = 10) 
    )
}




# ==============================================================================
# 3. 手动参数调整 (2x2 布局)
# ==============================================================================
print("正在绘制 2x2 面板...")

# --- 图 1: Hypoxic Microenvironment (新增) ---
# 缺氧评分通常较高，范围预估 -0.5 到 2.0
p_hypoxia <- plot_box_strict_manual(
  plot_data, "Hypoxia", "Hypoxic Microenvironment", 
  star_y = 1.6,       
  view_min = -1,    
  view_max = 1.8      
)

# --- 图 2: Glycolysis (保持) ---
p_glyco <- plot_box_strict_manual(
  plot_data, "Glycolysis", "Glycolysis Score", 
  star_y = 1.7,       
  view_min = -0.8,    
  view_max = 1.9      
)

# --- 图 3: Lipid Biosynthesis (保持) ---
p_lipid_syn <- plot_box_strict_manual(
  plot_data, "Lipid_Synthesis", "Lipid Biosynthesis Score",
  star_y = 0.65,      
  view_min = -0.5, 
  view_max = 0.75
)

# --- 图 4: Lipid Uptake (保持你调整好的参数) ---
p_lipid_uptake <- plot_box_strict_manual(
  plot_data, "Lipid_Uptake", "Lipid Uptake Score",
  star_y = 0.12,      # 之前调好的完美高度
  view_min = -0.15,    
  view_max = 0.2     
)

# ==============================================================================
# 4. 拼图保存 (2行2列)
# ==============================================================================
# 逻辑布局：
# 上排：环境与能量 (缺氧 + 糖酵解) -> 都是 High 高
# 下排：脂质代谢 (合成高 + 摄取低) -> 形成对比

pdf("Fig5_Metabolic_Panel_2x2_Final.pdf", width = 10, height = 9)

# 使用 patchwork 拼图
layout <- (p_hypoxia + p_glyco) / (p_lipid_syn + p_lipid_uptake)

print(layout)
dev.off()

print("✅ 2x2 代谢全景图已生成：Fig5_Metabolic_Panel_2x2_Final.pdf")
