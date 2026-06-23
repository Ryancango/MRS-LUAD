# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/CellChat_Analysis/cellchat_circleplot.R
# Role: Draw CellChat circle plots for selected signalling pathways.
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
library(RColorBrewer)

# 设置路径
work_dir <- "~/lung_cancer/scRNA/131907/CellChat_Analysis"
setwd(work_dir)
cellchat <- readRDS("cellchat_object_computed.rds")

# ==============================================================================
# 1. 样式设置 (顶刊级配色)
# ==============================================================================
# 提取所有细胞名
group.cell <- levels(cellchat@idents)

# 定义 NPG 风格配色
my_colors <- scales::hue_pal()(length(group.cell))
names(my_colors) <- group.cell

# 核心配色 (High=红, Low=蓝)
my_colors["Epi_High"] <- "#E64B35" 
my_colors["Epi_Low"]  <- "#4DBBD5" 
# 接收者配色
my_colors["Myeloid cells"] <- "#8491B4"
my_colors["T lymphocytes"] <- "#91D1C2"
my_colors["Fibroblasts"] <- "#F39B7F"
my_colors["Endothelial cells"] <- "#00A087"

# ==============================================================================
# 2. 绘制“合二为一”的圆圈图
# ==============================================================================
# 设定通路
target_pathway <- "MIF"

# 设定发送者：同时包含 High 和 Low
sources_use <- c("Epi_High", "Epi_Low")

# 设定接收者
targets_use <- c("Myeloid cells", "T lymphocytes", "Fibroblasts", "Endothelial cells")

# 检查通路是否存在
if(target_pathway %in% cellchat@netP$pathways) {
  
  print(paste("正在绘制", target_pathway, "的合并圆圈图..."))
  
  pdf(paste0("Fig_Merged_", target_pathway, "_Circle.pdf"), width = 8, height = 8)
  
  # 【核心函数】 netVisual_aggregate
  # 关键点：sources.use 里同时写了两个 Epi
  # 这样它们就会出现在同一个圆圈里，发出的线条颜色会自动对应细胞颜色
  
  netVisual_aggregate(cellchat, 
                      signaling = target_pathway, 
                      sources.use = sources_use,
                      targets.use = targets_use,
                      layout = "circle", 
                      color.use = my_colors,
                      
                      # 美化参数
                      edge.width.max = 15,    # 线条加粗，方便对比
                      vertex.weight = 20,     # 节点点变大
                      weight.scale = TRUE,    # 根据强度缩放线条
                      label.edge = FALSE,     # 不显示具体的数字，保持画面干净
                      pt.title = 16,          # 标题大小
                      vertex.label.cex = 1.2  # 细胞标签大小
  )
  
  # 手动加个标题
  title(paste0("Merged ", target_pathway, " Signaling: High (Red) vs Low (Blue)"))
  
  dev.off()
  print("✅ 合并圆圈图已生成！")
  
} else {
  print("错误：该通路不存在。")
}