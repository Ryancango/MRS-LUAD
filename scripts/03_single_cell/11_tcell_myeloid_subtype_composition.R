# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/T细胞髓样细胞亚群.R
# Role: Summarize T-cell and myeloid-cell subtype composition by sample-level epithelial-risk grouping.
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
library(dplyr)
library(ggsci)
library(patchwork)

print(paste0("[", Sys.time(), "] 正在加载数据..."))
seurat_obj <- readRDS("seurat_tlung_riskgroup.rds")

# ==============================================================================
# 1. 定义通用函数：重新降维与绘图
# ==============================================================================
# 这个函数负责：提取子集 -> 过滤垃圾 -> 重新跑PCA/tSNE -> 画美图
process_and_plot <- function(obj, lineage_name, title_text, junk_types) {
  
  print(paste("正在处理:", title_text, "..."))
  
  # 1. 提取特定谱系
  sub_obj <- subset(obj, subset = Cell_type == lineage_name)
  
  # 2. 数据清洗 (剔除 NA, Undetermined, Mixed 等垃圾亚群)
  # 这样图会非常干净
  sub_obj <- subset(sub_obj, subset = Cell_subtype %in% junk_types, invert = TRUE)
  sub_obj <- subset(sub_obj, subset = !is.na(Cell_subtype))
  
  # 3. 重新标准化与降维 (关键步骤！让亚群聚成漂亮的形状)
  # Standard Seurat Workflow
  sub_obj <- NormalizeData(sub_obj) %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA(verbose = FALSE) %>% 
    RunTSNE(dims = 1:15) # 重新计算 tSNE 坐标
  
  # 4. 准备配色 (NPG 风格扩展)
  n_groups <- length(unique(sub_obj$Cell_subtype))
  my_cols <- colorRampPalette(pal_npg("nrc")(10))(n_groups)
  
  # 5. 绘图
  p <- DimPlot(sub_obj, 
               reduction = "tsne", 
               group.by = "Cell_subtype", 
               cols = my_cols,
               pt.size = 0.5,     # 点稍微大一点，更有质感
               label = TRUE,      # 显示标签
               label.size = 4,    # 标签字号
               label.box = TRUE,  # 标签加白底框
               repel = TRUE) +    # 防重叠
    
    # 极简主题
    NoAxes() + 
    ggtitle(title_text) +
    
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right",
      legend.text = element_text(size = 10),
      aspect.ratio = 1 # 保持正方形
    )
  
  return(p)
}

# ==============================================================================
# 2. 执行绘图
# ==============================================================================

# --- A. 处理 T 细胞 ---
# 定义要剔除的“垃圾”亚群 (参考你之前的 Cleaned 图)
t_junk <- c("NA", "Undetermined", "CD8+/CD4+ Mixed Th", "NK") # NK通常单独画，或者你可以保留

p_t <- process_and_plot(seurat_obj, "T lymphocytes", "T Cell Subtypes Landscape", t_junk)

# --- B. 处理 Myeloid 细胞 ---
m_junk <- c("NA", "Undetermined")

p_m <- process_and_plot(seurat_obj, "Myeloid cells", "Myeloid Subtypes Landscape", m_junk)

# ==============================================================================
# 3. 拼图与保存
# ==============================================================================
print("正在拼图...")

# 1行2列
final_plot <- p_t + p_m + plot_layout(ncol = 2)

pdf("Fig7_Subtype_Clustering_Panel.pdf", width = 16, height = 7)
print(final_plot)
dev.off()

print("✅ 聚类亚群图已生成：Fig7_Subtype_Clustering_Panel.pdf")