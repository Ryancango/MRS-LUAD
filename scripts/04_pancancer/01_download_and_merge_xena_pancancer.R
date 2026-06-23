# Curated release script for: MRS-LUAD
# Source script: intersection_mechine/mechine/outside_data/xena/download.R
# Role: Download and merge TCGA Pan-Cancer expression and survival endpoints from UCSC Xena.
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
# 0. 设置你的工作目录和目标路径
# ==============================================================================
library(data.table)
library(dplyr)

# 这是你刚才提供的自定义目录
work_dir <- file.path(mrs_project_root(), "intersection_mechine/mechine/outside_data/xena/")

# 确保目录存在，不存在则创建
if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
}

# 目标文件路径
target_file <- file.path(work_dir, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")

# ==============================================================================
# 1. 清理坏文件并重新下载 (使用 wget)
# ==============================================================================
# 检查是否存在旧文件，如果存在且小于 300MB，则视为损坏并删除
if (file.exists(target_file)) {
  fsize <- file.info(target_file)$size / 1024 / 1024
  if (fsize < 300) {
    print(paste("发现损坏的旧文件 (", round(fsize, 2), "MB)，正在删除..."))
    file.remove(target_file)
  } else {
    print("检测到本地已有完整文件 (300MB+)，跳过下载。")
  }
}

# 如果文件不存在（被删除了或本来就没有），则开始下载
if (!file.exists(target_file)) {
  url <- "https://pancanatlas.xenahubs.net/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
  
  # 构造 wget 命令 (Linux系统最稳健的下载方式)
  # -O 指定输出文件名
  # -c 支持断点续传
  cmd <- paste0('wget -c -O "', target_file, '" "', url, '"')
  
  print("正在开始下载 (文件较大，约315MB，请耐心等待)...")
  system(cmd)
}

# ==============================================================================
# 2. 验证下载结果
# ==============================================================================
final_size <- file.info(target_file)$size / 1024 / 1024
print(paste("当前文件大小:", round(final_size, 2), "MB"))

if (final_size < 300) {
  stop("下载似乎又失败了（文件太小）。请检查网络，或尝试在终端手动运行 wget 命令。")
}

print("文件完整！正在读取数据...")

# ==============================================================================
# 3. 读取与转置 (Processing)
# ==============================================================================
rna_expr_data <- fread(target_file)

print(paste("原始数据维度:", paste(dim(rna_expr_data), collapse = " x ")))
# 此时应该是 20531 x 11070

print("正在转置矩阵 (行变样本，列变基因)...")
# 提取基因名
gene_names <- rna_expr_data[[1]]
# 转置数值矩阵
expr_mat_t <- t(as.matrix(rna_expr_data[, -1]))
# 赋列名
colnames(expr_mat_t) <- gene_names
# 转回数据框
expr_df <- as.data.frame(expr_mat_t)
# 添加 ID 列
expr_df$ID <- rownames(expr_mat_t)

# 清理内存
rm(rna_expr_data, expr_mat_t)
gc()

# ==============================================================================
# 4. 整理临床数据并合并
# ==============================================================================
# 假设 survival_data 还在你的环境中，如果不在，请先加载它
if(!exists("survival_data")) {
  # 如果你没有存 survival_data，这里尝试重新下载（因为这个文件很小，下载很快）
  print("未找到生存数据，正在重新获取...")
  library(UCSCXenaTools)
  surv_info <- XenaGenerate() %>%
    XenaFilter(filterDatasets = "Survival_SupplementalTable_S1_20171025_xena_sp") %>%
    XenaQuery() %>%
    XenaDownload()
  survival_data <- XenaPrepare(surv_info)
}

print("正在整理临床数据...")
surv_clean <- survival_data %>%
  dplyr::select(
    sample, 
    OS.time, OS, 
    DSS.time, DSS, 
    DFI.time, DFI, 
    PFI.time, PFI, 
    `cancer type abbreviation`
  ) %>%
  dplyr::rename(
    ID = sample, 
    CancerType = `cancer type abbreviation` # 改英文名
  ) %>%
  distinct(ID, .keep_all = TRUE)

print("正在合并...")
final_data <- merge(surv_clean, expr_df, by = "ID")

# ==============================================================================
# 5. 保存最终结果
# ==============================================================================
print(paste("最终数据维度:", paste(dim(final_data), collapse = " x ")))

# 保存路径
save_path <- file.path(work_dir, "TCGA_PANCAN_Merged_Final.rds")

print(paste("正在保存结果至:", save_path))
saveRDS(final_data, file = save_path)

print("全部完成！")