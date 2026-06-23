# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/CellChat_Analysis/小提琴堆叠图.R
# Role: Plot expression of selected ligand-receptor genes across cell groups.
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
library(patchwork)

# 设置路径 (请确认路径无误)
work_dir <- "~/lung_cancer/scRNA/131907/CellChat_Analysis"
setwd(work_dir)

print(paste0("[", Sys.time(), "] 正在加载数据..."))
cellchat <- readRDS("cellchat_object_computed.rds")

# ==============================================================================
# 1. 定义基因列表 (按 机制模块 排序)
# ==============================================================================
# 这里的逻辑是：先放 High Risk 特有的配体，再放对应的受体
# 这样从上往下看，逻辑非常顺畅

genes_to_plot <- c(
  # --- 1. 免疫抑制 (MIF Axis) ---
  "MIF",             # 配体：预期 Epi_High 高
  "CD74", "CXCR4",   # 受体：预期 Myeloid/T cells 高
  "CD44",            # 辅助受体
  
  # --- 2. 血管生成 (VEGF Axis) ---
  "VEGFB",           # 配体：预期 Epi_High 高
  "FLT1",            # 受体 (即 VEGFR1)：预期 Endothelial 高
  
  # --- 3. 基质重塑 (Endothelin Axis) ---
  "EDN1",            # 配体：预期 Epi_High 高
  "EDNRA",           # 受体：预期 Fibroblasts 高
  
  # --- 4. 免疫检查点 (Galectin Axis) ---
  "LGALS9",          # 配体：预期 Epi_High 高
  "HAVCR2",          # 受体 (TIM-3)：预期 T cells/Myeloid 高
  "P4HB"             # 受体
)

# 检查基因是否存在于数据矩阵中 (防止报错)
available_genes <- rownames(cellchat@data.signaling)
# 取交集
valid_genes <- intersect(genes_to_plot, available_genes)
# 保持我们定义的顺序
valid_genes <- genes_to_plot[genes_to_plot %in% valid_genes]

print("即将绘制以下基因:")
print(valid_genes)

# ==============================================================================
# 2. 定义统一配色 (NPG Style)
# ==============================================================================
group.cell <- levels(cellchat@idents)
my_colors <- scales::hue_pal()(length(group.cell))
names(my_colors) <- group.cell

# 核心配色锁定
my_colors["Epi_High"] <- "#E64B35" # 红色
my_colors["Epi_Low"]  <- "#4DBBD5" # 蓝色
# 靶细胞配色
if("Myeloid cells" %in% names(my_colors)) my_colors["Myeloid cells"] <- "#8491B4" # 紫色
if("T lymphocytes" %in% names(my_colors)) my_colors["T lymphocytes"] <- "#91D1C2" # 青色
if("Fibroblasts" %in% names(my_colors))   my_colors["Fibroblasts"] <- "#F39B7F"   # 橙色
if("Endothelial cells" %in% names(my_colors)) my_colors["Endothelial cells"] <- "#00A087" # 深绿

# ==============================================================================
# 3. 绘制堆叠小提琴图
# ==============================================================================
print("正在绘图...")

# 计算合适的高度：每个基因给 0.8 inch，加上标题空间
plot_height <- length(valid_genes) * 0.8 + 2

pdf("Fig6_Detail_Stacked_Violin_Mechanisms.pdf", width = 8, height = plot_height)

p <- plotGeneExpression(cellchat, 
                        features = valid_genes, 
                        group.by = "cellchat_group", 
                        color.use = my_colors,
                        type = "violin") +
  
  # --- 美化细节 ---
  theme(
    # X轴文字：倾斜45度，加粗，黑色
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
    
    # 基因名称(Y轴)：加粗，水平放置 (angle=0)
    strip.text.y = element_text(angle = 0, size = 10, face = "bold", hjust = 0, color = "black"),
    
    # 去掉条带背景，改用干净的白底
    strip.background = element_rect(fill = "white", color = NA), 
    
    # 去掉杂乱的背景颜色，只留虚线辅助
    panel.background = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  ) +
  
  # 添加标题
  labs(title = "Molecular Mechanism Validation",
       subtitle = "Expression of Key Ligands (High-Risk) and Receptors (TME)")

print(p)
dev.off()

print(paste0("✅ 成功！文件已保存至: ", work_dir, "/Fig6_Detail_Stacked_Violin_Mechanisms.pdf"))