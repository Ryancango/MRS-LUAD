# Curated release script for: MRS-LUAD
# Source script: wgcna分析/wgcna全基因.R
# Role: Run WGCNA on the batch-corrected GEO discovery matrix to identify tumour-associated modules.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}



# 清除工作区
rm(list = ls())


# 1. 加载必要的包

library(dplyr) # 加载dplyr包
library(WGCNA)
library(stringr)
library(ggplot2)
library(reshape2)
library(patchwork)

# 2. 数据读取与整合
expr_data2 <- read.csv("batch_corrected_expression_clean.csv", 
                       row.names = 1, check.names = FALSE)
dim(expr_data2)  # 显示基因数×样本数

expr_data <- t(expr_data2)

sample_names <- colnames(expr_data)
head(sample_names, 5)

group <- ifelse(str_detect(sample_names, "Control"), "Control", "Tumor")
table(group)

traits <- data.frame(sample = sample_names, group = group)
rownames(traits) <- sample_names

# 3. 检查基因和样本质量
datExpr <- t(expr_data)  # 转为样本×基因
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 4. 方差筛选前5000基因
gene_variances <- apply(datExpr, 2, var)
n_genes <- min(12000, ncol(datExpr))

top_genes <- names(sort(gene_variances, decreasing = TRUE)[1:n_genes])
datExpr_filtered <- datExpr[, top_genes]
dim(datExpr_filtered)
#datExpr_filtered <- t(datExpr_filtered)

# 5. WGCNA分析
# (1) 软阈值选择
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr_filtered, powerVector = powers, verbose = 5)
selected_power <- sft$powerEstimate
fit_indices <- data.frame(
  Power = powers,
  Signed.R.sq = -sign(sft$fitIndices[, "slope"]) * sft$fitIndices[, "SFT.R.sq"],
  mean.k = sft$fitIndices[, "mean.k."]
)

scale_independence_plot <- ggplot(fit_indices, aes(x = Power, y = Signed.R.sq)) +
  geom_text(aes(label = Power), color = "red", size = 4) +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = selected_power, color = "blue", linetype = "dotted", linewidth = 0.8) +
  labs(x = "Soft Threshold (power)", y = "Scale-free Topology Model Fit (signed R²)") +
  theme_bw()

mean_connectivity_plot <- ggplot(fit_indices, aes(x = Power, y = mean.k)) +
  geom_text(aes(label = Power), color = "blue", size = 4) +
  geom_vline(xintercept = selected_power, color = "blue", linetype = "dotted", linewidth = 0.8) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Soft Threshold (power)", y = "Mean Connectivity (log scale)") +
  theme_bw()

if (is.na(selected_power)) {
  warning("Automatic power selection failed, please choose manually!")
  selected_power <- 6
}

combined_plots <- scale_independence_plot + mean_connectivity_plot +
  plot_annotation(title = "WGCNA Network Topology Analysis")
print(combined_plots)

# (2) 构建基因网络
#softPower <- selected_power

# 强制使用 WGCNA 包的 cor 函数
cor <- WGCNA::cor
net <- blockwiseModules(
  datExpr_filtered,
  power = 4,
  maxBlockSize = 40000,
  TOMType = "unsigned",
  minModuleSize = 50,
  deepSplit = 2,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  corType = "pearson",
  verbose = 3
)
# 恢复 R 默认的 cor 函数
cor <- stats::cor

module_colors <- labels2colors(net$colors)
plotDendroAndColors(
  net$dendrograms[[1]],
  module_colors,
  "Module Colors",
  dendroLabels = FALSE,
  hang = 0.03)

# 6. 表型关联分析
trait_data <- data.frame(
  Control = as.numeric(grepl("^Control", traits$group)),
  Tumor = as.numeric(grepl("^Tumor", traits$group)),
  row.names = rownames(traits)
)

MEs <- moduleEigengenes(datExpr_filtered, module_colors)$eigengenes
moduleTraitCor <- cor(MEs, trait_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr_filtered))


# 7. 可视化结果
heatmap_data <- expand.grid(
  Module = colnames(MEs),
  Trait = colnames(trait_data)
) %>%
  mutate(
    Correlation = as.vector(moduleTraitCor),
    Pvalue = as.vector(moduleTraitPvalue)
  ) %>%
  mutate(Significance = case_when(
    Pvalue < 0.001 ~ "***",
    Pvalue < 0.01 ~ "**",
    Pvalue < 0.05 ~ "*",
    TRUE ~ ""
  ))

ggplot(heatmap_data, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(Correlation, 2), "\n", Significance)), 
            color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Module-Trait Relationships", x = "Experimental Traits", y = "WGCNA Modules")



# 计算模块成员度（MM）
MM <- cor(datExpr_filtered, MEs, use = "p")

# 计算基因显著性（GS），以Tumor为例作为表型变量
GS <- cor(datExpr_filtered, trait_data$Tumor, use = "p")

# 获取模块列表
modules <- unique(module_colors)

# 为每个模块绘制散点图
for (module in modules) {
  # 获取模块中的基因
  module_genes <- which(module_colors == module)
  
  # 获取该模块的MM和GS
  MM_module <- MM[module_genes, paste0("ME", module)]
  GS_module <- GS[module_genes]
  
  # 计算MM和GS的相关性及其P值
  cor_MM_GS <- cor(MM_module, GS_module, use = "p")
  p_value <- corPvalueStudent(cor_MM_GS, length(module_genes))
  
  # 准备绘图数据
  plot_data <- data.frame(MM = MM_module, GS = GS_module)
  
  # 绘制美化后的散点图
  p <- ggplot(plot_data, aes(x = MM, y = GS)) +
    geom_point(color = "steelblue", alpha = 0.7, size = 2) +  # 散点样式
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 回归线
    labs(title = paste("Module", module, "MM vs GS Scatter Plot"),  # 标题
         x = "Module Membership (MM)",  # X轴标签
         y = "Gene Significance (GS)") +  # Y轴标签
    annotate("text", 
             x = min(MM_module) + 0.05 * diff(range(MM_module)), 
             y = max(GS_module) - 0.05 * diff(range(GS_module)),
             label = paste("Correlation:", round(cor_MM_GS, 2), 
                           "\nP-value:", format.pval(p_value, digits = 2)),
             hjust = 0, vjust = 1, size = 4, color = "black") +  # 注释
    theme_bw() +  # 黑白主题
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题样式
      axis.title = element_text(size = 12),  # 轴标题字体
      axis.text = element_text(size = 10),  # 轴刻度字体
      panel.grid.major = element_blank(),  # 去掉主网格线
      panel.grid.minor = element_blank()  # 去掉次网格线
    )
  
  # 显示图表
  print(p)
}

