# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/CellChat_Analysis/cellchatplot1.R
# Role: Plot CellChat interaction-count and interaction-strength summaries.
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
library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)

# 设置路径
work_dir <- "~/lung_cancer/scRNA/131907/CellChat_Analysis"
setwd(work_dir)

print(paste0("[", Sys.time(), "] 正在加载数据..."))
cellchat <- readRDS("cellchat_object_computed.rds")

# 1. 提取通讯矩阵
mat <- cellchat@net$weight

# 2. 定义颜色方案 (保持 NPG 风格)
cell_names <- rownames(mat)
col_cells <- scales::hue_pal()(length(cell_names))
names(col_cells) <- cell_names
col_cells["Epi_High"] <- "#E64B35" # 红
col_cells["Epi_Low"]  <- "#4DBBD5" # 蓝
other_cells <- setdiff(cell_names, c("Epi_High", "Epi_Low"))
col_cells[other_cells] <- c("#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#7E6148", "#DC0000")[1:length(other_cells)]

# ==============================================================================
# 1. 顶刊级配色 (烈焰金方案)
# ==============================================================================
# 白色底 -> 淡金 -> 亮红 -> 深紫
col_fun <- colorRamp2(c(0, max(mat)/4, max(mat)/2, max(mat)), 
                      c("#FFFFFF", "#FFEDA0", "#E03E3E", "#4A0C44"))

# ==============================================================================
# 2. 构建注释条 (精细化字体与标签)
# ==============================================================================
row_sums <- rowSums(mat)
col_sums <- colSums(mat)

# --- 顶部注释 (Total Incoming) ---
top_anno <- HeatmapAnnotation(
  # 1. 柱状图
  "Total Incoming" = anno_barplot(col_sums, 
                                  gp = gpar(fill = col_cells[names(col_sums)], col = NA), 
                                  border = FALSE, 
                                  height = unit(2, "cm"),
                                  # 【字体修正】加粗轴标签，字号10
                                  axis_param = list(gp = gpar(fontsize = 10, fontface = "bold"))),
  
  # 2. 颜色块
  "Cell Type" = anno_simple(colnames(mat), col = col_cells),
  
  col = list("Cell Type" = col_cells),
  
  # 【核心修改】只显示第一个(Total Incoming)，隐藏第二个(Cell Type)
  show_annotation_name = c(TRUE, FALSE), 
  
  # 名字放在左边
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  simple_anno_size = unit(0.4, "cm")
)

# --- 右侧注释 (Total Outgoing) ---
right_anno <- rowAnnotation(
  # 1. 颜色块
  "Cell Type" = anno_simple(rownames(mat), col = col_cells),
  
  # 2. 柱状图
  "Total Outgoing" = anno_barplot(row_sums, 
                                  gp = gpar(fill = col_cells[names(row_sums)], col = NA), 
                                  border = FALSE, 
                                  width = unit(2, "cm"),
                                  # 【字体修正】加粗轴标签
                                  axis_param = list(gp = gpar(fontsize = 10, fontface = "bold"))),
  
  col = list("Cell Type" = col_cells),
  
  # 【核心修改】只显示第二个(Total Outgoing)，隐藏第一个(Cell Type)
  show_annotation_name = c(FALSE, TRUE),
  
  # 名字放在底部
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

# ==============================================================================
# 3. 绘制最终图
# ==============================================================================
pdf("Fig6_Overview_Heatmap_Final_Clean_v2.pdf", width = 10, height = 8.5)

Heatmap(mat,
        name = "Interaction\nStrength", 
        col = col_fun,
        
        # --- 边框美化 ---
        rect_gp = gpar(col = "white", lwd = 2), # 瓷砖质感
        
        # --- 聚类设置 ---
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        
        # --- 轴标签美化 (统一加粗) ---
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "sans"),
        
        column_names_side = "bottom",
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "sans"),
        
        # --- 绑定注释 ---
        top_annotation = top_anno,
        right_annotation = right_anno,
        
        # --- 标题 ---
        column_title = "Receiver (Target)",
        row_title = "Sender (Source)",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        row_title_gp = gpar(fontsize = 14, fontface = "bold")
)

dev.off()
print("✅ 完美版热图已生成：Fig6_Overview_Heatmap_Final_Clean_v2.pdf")