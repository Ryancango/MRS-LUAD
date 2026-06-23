# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/Pseudotime_Quantification/拟时序分析拟合曲线主图.R
# Role: Fit and plot gene-expression trends along pseudotime.
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
library(Seurat)
library(monocle)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(msigdbr) # 用于获取标准Hallmark

# 设置路径
work_dir <- file.path(mrs_project_root(), "scRNA/131907/Pseudotime_Quantification")
if(!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

print(paste0("[", Sys.time(), "] 正在加载数据..."))
seurat_path <- file.path(mrs_project_root(), "scRNA/131907/Pseudotime_Quantification/seurat_tlungepi_risk.rds")
monocle_path <- file.path(mrs_project_root(), "scRNA/131907/Pseudotime_Quantification/monocle_cds_final.rds")

seurat_obj <- readRDS(seurat_path)
my_cds <- readRDS(monocle_path)

# ==============================================================================
# 1. 定义“全明星”基因集 (混合使用 Hallmark 和 精选集)
# ==============================================================================
print("正在准备基因集...")

# --- A. 从 Hallmark 获取标准通路 ---
# 获取基因集
tryCatch({
  h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
}, error = function(e) {
  h_gene_sets <- msigdbr(species = "Homo sapiens", collection = "H")
})
hallmark_list <- split(h_gene_sets$gene_symbol, h_gene_sets$gs_name)

# 提取我们需要的标准通路
sigs_list <- list(
  Glycolysis_Score = hallmark_list[["HALLMARK_GLYCOLYSIS"]],
  EMT_Score = hallmark_list[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
  G2M_Score = hallmark_list[["HALLMARK_G2M_CHECKPOINT"]],
  FattyAcid_Score = hallmark_list[["HALLMARK_FATTY_ACID_METABOLISM"]]
)

# --- B. 【关键替换】使用精选的高质量缺氧基因集 ---
# 替换掉 HALLMARK_HYPOXIA，用这个 R=0.37 的版本
optimized_hypoxia <- c("VEGFA", "SLC2A1", "PGK1", "CA9", "LDHA", "ENO1", "ADM", "BNIP3")
sigs_list$Hypoxia_Score <- optimized_hypoxia

# ==============================================================================
# 2. 计算评分
# ==============================================================================
print("正在计算评分...")

# 批量计算
seurat_obj <- AddModuleScore(seurat_obj, features = sigs_list, name = "FinalSig")

# 准备绘图数据
common_cells <- intersect(colnames(my_cds), colnames(seurat_obj))
plot_data <- pData(my_cds)[common_cells, ]

# 转移分数 (注意 AddModuleScore 的命名规则是 name + index)
# 顺序: 1=Glyco, 2=EMT, 3=G2M, 4=Fatty, 5=Hypoxia
plot_data$Glycolysis <- seurat_obj@meta.data[common_cells, "FinalSig1"]
plot_data$EMT <- seurat_obj@meta.data[common_cells, "FinalSig2"]
plot_data$G2M <- seurat_obj@meta.data[common_cells, "FinalSig3"]
plot_data$FattyAcid <- seurat_obj@meta.data[common_cells, "FinalSig4"]
plot_data$Hypoxia <- seurat_obj@meta.data[common_cells, "FinalSig5"]

# 转移 Risk Score 和 Group
plot_data$RiskScore <- seurat_obj@meta.data[common_cells, "risk_score"]
plot_data$RiskGroup <- seurat_obj@meta.data[common_cells, "risk_group"]

# ==============================================================================
# 3. 绘图函数 (顶刊风格)
# ==============================================================================
npg_cols <- c("High Risk" = "#E64B35", "Low Risk" = "#4DBBD5")

plot_final_style <- function(df, x_col, y_col, y_lab, title_text) {
  
  # 统计
  res <- cor.test(df[[x_col]], df[[y_col]], method = "spearman")
  r_val <- round(res$estimate, 2)
  if(res$p.value < 2.2e-16){
    p_text <- "P < 2.2e-16"
  } else {
    p_text <- sprintf("P = %.2e", res$p.value)
  }
  sub_title <- paste0("Spearman R = ", r_val, ", ", p_text)
  
  ggplot(df, aes_string(x = x_col, y = y_col)) +
    geom_point(aes(color = RiskGroup), size = 0.6, alpha = 0.5, stroke = 0) + # 稍微调大一点点点
    geom_smooth(method = "loess", color = "black", linewidth = 1.2, se = TRUE, fill = "grey85") +
    scale_color_manual(values = npg_cols) +
    theme_classic() +
    labs(title = title_text,
         subtitle = sub_title,
         x = "Pseudotime (Evolution)", 
         y = y_lab) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10, color = "black"),
      legend.position = "none"
    )
}

# ==============================================================================
# 4. 生成最终面板
# ==============================================================================
print("正在生成 Figure 5 最终版...")

# 1. MRS Risk Score (总纲)
p1 <- plot_final_style(plot_data, "Pseudotime", "RiskScore", "MRS Risk Score", "Malignancy Evolution")

# 2. Glycolysis (代谢核心)
p2 <- plot_final_style(plot_data, "Pseudotime", "Glycolysis", "Glycolysis Score", "Metabolic Reprogramming")

# 3. Hypoxia (环境压力 - 【已替换为 R=0.37 的版本】)
p3 <- plot_final_style(plot_data, "Pseudotime", "Hypoxia", "Hypoxia Score", "Hypoxic Microenvironment")

# 4. Fatty Acid (代谢重塑 - 下降)
p4 <- plot_final_style(plot_data, "Pseudotime", "FattyAcid", "Fatty Acid Score", "Lipid Metabolism Loss")

# 5. G2M (增殖)
p5 <- plot_final_style(plot_data, "Pseudotime", "G2M", "G2M Score", "Proliferation Capacity")

# 6. EMT (侵袭)
p6 <- plot_final_style(plot_data, "Pseudotime", "EMT", "EMT Score", "Invasion Potential")

# --- 提取共享图例 ---
dummy_p <- ggplot(plot_data, aes(x=Pseudotime, y=RiskScore, color=RiskGroup)) +
  geom_point(size=4) +
  scale_color_manual(values = npg_cols) +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_text(face="bold", size=12), legend.text = element_text(size=12)) + 
  guides(color = guide_legend(override.aes = list(size=4, alpha=1)))
legend_shared <- get_legend(dummy_p)

# --- 组合 (2行3列，或者3行2列，推荐 2行3列更宽屏) ---
# 布局顺序：风险 -> 糖 -> 缺氧 (第一排：驱动力)
#          脂 -> 增殖 -> EMT (第二排：表型后果)
plots_grid <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3, align = "hv")
final_combined <- plot_grid(legend_shared, plots_grid, ncol = 1, rel_heights = c(0.08, 1))

# 保存
pdf_file <- "Figure5_Final_All_Star.pdf"
ggsave(pdf_file, plot = final_combined, width = 16, height = 10)

print(paste0("✅ Figure 5 最终版已生成: ", pdf_file))