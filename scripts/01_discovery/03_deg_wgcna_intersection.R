# Curated release script for: MRS-LUAD
# Source script: intersection_mechine/DEG和module取交集.R
# Role: Intersect recurrent differentially expressed genes with the tumour-associated WGCNA module.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}

# 加载必要的库
library(readr)
library(dplyr)

# 读取CSV文件
df1 <- read_csv('common_up_regulated_genes.csv', col_names = TRUE)
df2 <- read_csv('module_blue_genes.csv', col_names = TRUE)

# 提取基因symbol列（第一列），去除空值，转换为字符串并去除首尾空格
gene_set1 <- df1[[1]] %>% na.omit() %>% as.character() %>% trimws()
gene_set2 <- df2[[1]] %>% na.omit() %>% as.character() %>% trimws()

# 取两个基因列表的交集
intersection <- intersect(gene_set1, gene_set2)

# 将交集转换为数据框
intersection_df <- data.frame(gene_symbol = intersection)

# 保存交集到新的CSV文件
write_csv(intersection_df, 'blue&up_genes.csv')
