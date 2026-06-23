# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/T细胞耗竭细胞比例图.R
# Role: Quantify exhausted T-cell proportions across high- and low-risk samples.
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
# 0. 环境与数据
# ==============================================================================
rm(list = ls())
setwd("~/lung_cancer/scRNA/131907")
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(patchwork)

print(paste0("[", Sys.time(), "] 正在加载数据..."))
seurat_obj <- readRDS("seurat_tlung_riskgroup.rds")

# ==============================================================================
# 1. 恢复样本风险分组
# ==============================================================================
meta_epi <- subset(seurat_obj, subset = Cell_type == "Epithelial cells")@meta.data
sample_stats <- meta_epi %>%
  group_by(Sample) %>%
  summarise(n_high = sum(risk_group == "High Risk"), total = n(), ratio = n_high/total)
threshold <- median(sample_stats$ratio)
sample_stats$Sample_Group <- ifelse(sample_stats$ratio > threshold, "High Risk", "Low Risk")
risk_map <- setNames(sample_stats$Sample_Group, sample_stats$Sample)
seurat_obj@meta.data$Sample_Group <- unname(risk_map[as.character(seurat_obj@meta.data$Sample)])
seurat_obj <- subset(seurat_obj, subset = !is.na(Sample_Group))

# ==============================================================================
# 2. 定义“黑名单” (清洗垃圾数据)
# ==============================================================================
junk_labels <- c("NA", "Undetermined", "CD8+/CD4+ Mixed Th", "NK")

# ==============================================================================
# 3. 绘制 T 细胞纯净版比例图 (修复变量名)
# ==============================================================================
print("正在绘制 T 细胞图...")
t_obj <- subset(seurat_obj, subset = Cell_type == "T lymphocytes")
df_t <- t_obj@meta.data

# 清洗
df_t_clean <- df_t %>% 
  filter(!Cell_subtype %in% junk_labels) %>% 
  filter(!is.na(Cell_subtype))

# 统计
data_t_summary <- df_t_clean %>%
  group_by(Sample_Group, Cell_subtype) %>%
  summarise(Count = n()) %>%
  mutate(Freq = Count / sum(Count)) %>%
  ungroup()

# 排序 (修复点：这里使用 data_t_summary)
order_t <- data_t_summary %>% 
  filter(Sample_Group == "High Risk") %>% 
  arrange(desc(Freq)) %>% 
  pull(Cell_subtype)

# 重新设置因子水平
data_t_summary$Cell_subtype <- factor(data_t_summary$Cell_subtype, levels = order_t)

# 绘图函数
plot_clean_bar <- function(data, title_text) {
  n_cols <- length(unique(data$Cell_subtype))
  # 动态颜色
  my_pal <- if(n_cols <= 10) pal_npg("nrc")(n_cols) else colorRampPalette(pal_npg("nrc")(10))(n_cols)
  
  ggplot(data, aes(x = Sample_Group, y = Freq, fill = Cell_subtype)) +
    geom_bar(stat = "identity", position = "fill", width = 0.6, color = "black", size = 0.2) +
    scale_fill_manual(values = my_pal) +
    theme_classic() +
    labs(title = title_text, y = "Fraction within Lineage", x = "") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      legend.title = element_blank()
    )
}

p_t <- plot_clean_bar(data_t_summary, "T Cell Subtypes (Cleaned)")

# ==============================================================================
# 4. 绘制 Myeloid 纯净版比例图
# ==============================================================================
print("正在绘制 Myeloid 图...")
m_obj <- subset(seurat_obj, subset = Cell_type == "Myeloid cells")
df_m <- m_obj@meta.data

# 清洗
df_m_clean <- df_m %>% 
  filter(!Cell_subtype %in% c("NA", "Undetermined")) %>%
  filter(!is.na(Cell_subtype))

# 统计
data_m_summary <- df_m_clean %>%
  group_by(Sample_Group, Cell_subtype) %>%
  summarise(Count = n()) %>%
  mutate(Freq = Count / sum(Count)) %>%
  ungroup()

# 排序 (使用 data_m_summary)
order_m <- data_m_summary %>% 
  filter(Sample_Group == "High Risk") %>% 
  arrange(desc(Freq)) %>% 
  pull(Cell_subtype)

data_m_summary$Cell_subtype <- factor(data_m_summary$Cell_subtype, levels = order_m)

p_m <- plot_clean_bar(data_m_summary, "Myeloid Subtypes (Cleaned)")

# ==============================================================================
# 5. 保存
# ==============================================================================
pdf("Fig7_Composition_Cleaned_Fixed.pdf", width = 12, height = 6)
print(p_t + p_m)
dev.off()

print("✅ 修复完成！图已生成：Fig7_Composition_Cleaned_Fixed.pdf")

