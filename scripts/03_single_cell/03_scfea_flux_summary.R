# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/scfea/scFEA.R
# Role: Summarize scFEA-inferred metabolic flux modules across cell types and risk groups.
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
# 0. 环境与包加载
# ==============================================================================
rm(list = ls())
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 设置工作目录
work_dir <- "~/lung_cancer/scRNA/131907/scFEA_Analysis"
if(!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

# ==============================================================================
# 1. 确认并配置输入文件路径 (!!!请务必检查这里!!!)
# ==============================================================================
# 请将这里的路径替换为你服务器上真实的绝对路径
flux_file <- file.path(scfea_project_root(), "scFEA_output/predicted_flux.csv") 
seurat_path <- file.path(mrs_project_root(), "scRNA/131907/seurat_tlung_riskgroup.rds")

if(!file.exists(flux_file)) stop(paste("错误: 找不到 Flux 文件", flux_file))
if(!file.exists(seurat_path)) stop(paste("错误: 找不到 Seurat 文件", seurat_path))

# ==============================================================================
# 2. 读取数据与匹配
# ==============================================================================
print(paste0("[", Sys.time(), "] 正在加载数据..."))

# 1. 读取 Flux 矩阵
flux_mat <- read.csv(flux_file, row.names = 1, check.names = FALSE)

# 2. 读取 Seurat 对象并提取最新的 Metadata
seurat_obj <- readRDS(seurat_path)
meta_data <- seurat_obj@meta.data

print("正在对齐细胞...")

# 3. 处理可能存在的细胞名不一致 (比如 "-" 变成了 ".")
if(length(intersect(rownames(flux_mat), rownames(meta_data))) == 0) {
  print("检测到行名格式不一致，尝试替换分隔符...")
  rownames(flux_mat) <- gsub("\\.", "-", rownames(flux_mat))
}

common_cells <- intersect(rownames(flux_mat), rownames(meta_data))
print(paste("最终匹配到的细胞数量:", length(common_cells)))

if(length(common_cells) < 100) {
  stop("错误: 匹配到的细胞太少！请手动检查 flux.csv 的行名和 Seurat 的 colnames 格式差异。")
}

# 取交集并合并
flux_mat <- flux_mat[common_cells, ]
meta_data <- meta_data[common_cells, ]
merged_df <- cbind(meta_data, flux_mat)

# ==============================================================================
# 3. 聚合计算与矩阵构建
# ==============================================================================
print("正在计算各细胞类型的代谢均值...")

# 确保 risk_group 是一个有序因子
merged_df$risk_group <- factor(merged_df$risk_group, levels = c("High Risk", "Low Risk"))

# 【兼容性修复】：确保 Cell_type 列名正确
if(!"Cell_type" %in% colnames(merged_df)){
  possible <- grep("cell.?type|ident", colnames(merged_df), ignore.case=T, value=T)
  if(length(possible) > 0) colnames(merged_df)[colnames(merged_df) == possible[1]] <- "Cell_type"
}

# 按 Cell_type 和 risk_group 分组，计算所有 M_ 开头的代谢模块的平均值
avg_df <- merged_df %>%
  group_by(Cell_type, risk_group) %>%
  summarise(across(starts_with("M_"), mean, .names = "{col}"), .groups = "drop")

# 构建热图专用的矩阵: 行是 M_xx，列是 "CellType|RiskGroup"
sample_ids <- paste0(avg_df$Cell_type, "|", avg_df$risk_group)
plot_mat <- as.matrix(avg_df %>% select(starts_with("M_")))
rownames(plot_mat) <- sample_ids
plot_mat <- t(plot_mat) 

# Z-score 标准化 (让红蓝差异更明显)
plot_mat_scaled <- t(scale(t(plot_mat)))
plot_mat_scaled[is.na(plot_mat_scaled)] <- 0 

# ==============================================================================
# 4. 准备 scFEA M168 模块注释表
# ==============================================================================
get_m168_annotation <- function() {
  anno_list <- list(
    "Glycolysis" = 1:11, "TCA Cycle" = 12:25, "Pentose Phosphate" = 26:33,
    "Serine" = 34:39, "Glycine" = 40:42, "Cysteine" = 43:47, "Threonine" = 48:52,
    "Lysine" = 53:57, "Aspartate/Asparagine" = 58:65, "Alanine" = 66:68,
    "BCAA" = 69:79, "Methionine" = 80:85, "Tyrosine" = 86:95,
    "Tryptophan" = 96:104, "Phenylalanine" = 105:109, "Glutamate/Gln" = 110:117,
    "Arginine/Proline" = 118:126, "Histidine" = 127:130, "Urea Cycle" = 131:135,
    "Fatty Acid" = 136:148, "Pyruvate" = 149:150, "Purine/Pyrimidine" = 151:160,
    "Transport" = 161:168
  )
  do.call(rbind, lapply(names(anno_list), function(n) {
    data.frame(Module_ID = paste0("M_", anno_list[[n]]), Super_Module = n)
  }))
}
module_anno <- get_m168_annotation()

# 合并注释与矩阵
df_temp <- as.data.frame(plot_mat_scaled)
df_temp$Module_ID <- rownames(df_temp)
df_merged <- merge(df_temp, module_anno, by = "Module_ID", all.x = TRUE)
df_merged <- df_merged[!is.na(df_merged$Super_Module), ] 

# 筛选变异度最高的前 60 个模块 (画全部 168 个图会太挤)
data_cols <- setdiff(colnames(df_merged), c("Module_ID", "Super_Module"))
df_merged$Var <- apply(df_merged[, data_cols], 1, var)

df_top <- df_merged %>% 
  arrange(desc(Var)) %>% 
  head(60) %>%
  arrange(Super_Module) # 这一步很关键，让同类代谢聚在一起

final_mat <- as.matrix(df_top[, data_cols])
rownames(final_mat) <- df_top$Module_ID
row_super_modules <- df_top$Super_Module

# ==============================================================================
# 5. 绘制并保存热图
# ==============================================================================
print("正在绘制高清热图...")

# 解析列名
col_groups <- colnames(final_mat)
cell_types <- sapply(strsplit(col_groups, "\\|"), `[`, 1)
risks <- sapply(strsplit(col_groups, "\\|"), `[`, 2)

# 配色
risk_colors <- c("High Risk" = "#E64B35", "Low Risk" = "#4DBBD5") # NPG
cts <- unique(cell_types)
ct_colors <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(cts)), cts)
sms <- unique(row_super_modules)
sm_colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(sms)), sms)

# 顶部注释条
top_anno <- HeatmapAnnotation(
  Risk = risks,
  CellType = cell_types,
  col = list(Risk = risk_colors, CellType = ct_colors),
  show_legend = TRUE,
  annotation_name_side = "left"
)

# 左侧大类注释条
left_anno <- rowAnnotation(
  Metabolism = row_super_modules,
  col = list(Metabolism = sm_colors),
  show_annotation_name = FALSE, 
  width = unit(4, "mm")
)

# 绘图输出
pdf_name <- "scFEA_Flux_Heatmap_Final.pdf"
pdf(pdf_name, width = 16, height = 12)

Heatmap(
  final_mat,
  name = "Flux Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#FFFFFF", "#B2182B")), 
  
  # 分列 (High 一半，Low 一半)
  column_split = factor(risks, levels = c("High Risk", "Low Risk")), 
  column_gap = unit(4, "mm"),
  
  # 分行 (按代谢种类)
  row_split = row_super_modules,
  row_gap = unit(1, "mm"),
  
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_rot = 45,
  column_labels = cell_types,
  
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  
  row_title_rot = 0, 
  row_title_side = "left",
  row_title_gp = gpar(fontface = "bold", fontsize = 10),
  
  top_annotation = top_anno,
  left_annotation = left_anno,
  
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5)
)

dev.off()
print(paste0("✅ 恢复成功！热图已保存至: ", file.path(getwd(), pdf_name)))