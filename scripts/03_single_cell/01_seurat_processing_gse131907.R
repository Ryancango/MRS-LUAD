# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/seurat.R
# Role: Process the GSE131907 single-cell data using the original annotation framework.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}

# 加载必要包
library(Seurat)
library(readxl)
library(dplyr)   # 用于数据处理

# 运行前释放内存
gc(full = TRUE)

# 基因列表
genes <- c("RRAS", "GPM6B", "FMO3", "CHRDL1", "ADRB2", "LPL", "ANGPT1", "SPOCK2", "P2RY14", "ADH1B")


# 运行SCTransform（完整数据集）
seurat_lung <- SCTransform(
  seurat_lung,
  method = "glmGamPoi",
  ncells = min(10000, ncol(seurat_lung)),  # 动态设置采样细胞数
  variable.features.n = 5000,              # 增加HVG数量
  vars.to.regress = NULL,          # 回归线粒体基因
  conserve.memory = TRUE,
  verbose = TRUE
)

# PCA降维
seurat_lung <- RunPCA(seurat_lung, assay = "SCT")

# 聚类分析
seurat_lung <- FindNeighbors(seurat_lung, dims = 1:30)
seurat_lung <- FindClusters(seurat_lung, resolution = 0.5)

# UMAP可视化
seurat_lung <- RunUMAP(seurat_lung, dims = 1:30)
DimPlot(seurat_lung, reduction = "umap")
DimPlot(seurat_lung, reduction = "tsne")




FeaturePlot(seurat_L, reduction = "tsne", features ="PHLDA2", split.by = "Sample_Origin", order = T)
FeaturePlot(seurat_L, reduction = "tsne", features ="PHLDA2", split.by = "Cell_type", order = T)

FeaturePlot(seurat_lung, reduction = "tsne", features ="RRAS", order = T)





# 绘制UMAP聚类图（按cell_type，无中文）
DimPlot(seurat_40k, reduction = "tsne", group.by = "Cell_type", label = TRUE, pt.size = 0.5)

DimPlot(seurat_LUNG, reduction = "tsne", group.by = "Sample_Origin", label = TRUE, pt.size = 0.5)


DimPlot(seu, reduction = "tsne", group.by = "risk_group", label = TRUE, pt.size = 0.1)



# VlnPlot绘图（ANOS1在cell_type中表达，按Type分类，无中文）
VlnPlot(seurat_lung, features = "ANOS1", group.by = "Cell_type", split.by = "Type")

# 基因列表
#genes <- c("EXOSC7", "PTPRB", "MAFF", "MT1M", "C5AR1", "PTRF", "TMEM47", "OLFML1", "ANOS1", "RGCC", "LRRC32", "ADGRE5")
genes <- c("RRAS", "GPM6B", "FMO3", "CHRDL1", "ADRB2", "LPL", "ANGPT1", "SPOCK2", "P2RY14", "ADH1B")
# 绘制堆积小提琴图（多基因堆积，按cell_type分组，无中文）
VlnPlot(seurat_40k, features = genes, group.by = "Cell_type", pt.size = 1, stack = TRUE)

