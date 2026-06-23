# Curated release script for: MRS-LUAD
# Source script: rawdata/预处理+PCA.R
# Role: Preprocess five GEO discovery cohorts, remove batch effects with ComBat, and inspect PCA outliers.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}

# 加载必要的包
library(sva)
library(tibble)
library(ggrepel)


# 文件列表
count_files <- c(
  "GSE7670_gene_symbol_expression_matrix.csv", 
  "GSE18842_gene_symbol_expression_matrix.csv",
  "GSE19804_gene_symbol_expression_matrix.csv", 
  "GSE31547_gene_symbol_expression_matrix.csv",
  #"GSE31552_gene_symbol_expression_matrix.csv", 
  "GSE33532_gene_symbol_expression_matrix.csv"
)

# 1. 读取并合并数据
merged_data <- NULL
batch_info <- list()
gene_intersect <- NULL

for (i in seq_along(count_files)) {
  # 读取CSV文件
  data <- read.csv(count_files[i], row.names = 1, check.names = FALSE)
  
  # 提取数据集ID作为批次前缀
  dataset_id <- gsub("_.*", "", basename(count_files[i]))
  
  # 重命名样本列名 (避免重复)
  colnames(data) <- paste0(dataset_id, "_", colnames(data))
  
  # 更新基因交集
  if (is.null(gene_intersect)) {
    gene_intersect <- rownames(data)
  } else {
    gene_intersect <- intersect(gene_intersect, rownames(data))
  }
  
  # 存储数据及批次信息
  batch_info[[i]] <- data
}

# 2. 筛选共同基因并合并
combined <- do.call(cbind, lapply(batch_info, function(df) {
  df[gene_intersect, , drop = FALSE]
}))

# 3. 转置矩阵 (样本为行，基因为列)
transposed_data <- t(combined)

# 4. 创建批次向量
batch_vector <- rep(seq_along(count_files), times = sapply(batch_info, ncol))

# 5. 使用ComBat进行去批次处理
# 注意：确保数据已标准化（通常需要log2转换）
# 这里假设数据已标准化，若未标准化需先进行log2(expr + 1)转换

# 创建模型矩阵（无协变量时使用空模型）
mod <- model.matrix(~1, data = as.data.frame(transposed_data))


transposed_data <- log2(transposed_data + 1)
# 执行ComBat去批次
cleaned_data <- ComBat(
  dat = t(transposed_data),  # 需要基因×样本的矩阵
  batch = batch_vector,
  mod = mod,
  par.prior = TRUE,   # 使用参数先验
  prior.plots = FALSE  # 不生成先验图
)

# 将结果转置回样本×基因格式
final_expression <- t(cleaned_data)

# 保存完整去批次矩阵（含离群值）
#write.csv(final_expression, "batch_corrected_expression_all_samples_without31552.csv")

#final_expression <- read.csv("batch_corrected_expression_all_samples_without31552.csv")

# 查看最终矩阵维度
dim(final_expression)

pca_res <- prcomp(final_expression, scale. = TRUE)

plot(pca_res$x[, 1:2], col = batch_vector, pch = 19)



# 7. PCA分析 (带标签)
pca_result <- prcomp(final_expression, scale. = TRUE)

# 创建PCA结果数据框
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Sample <- rownames(pca_df)
pca_df$Dataset <- sapply(strsplit(pca_df$Sample, "_"), `[`, 1)
pca_df$Batch <- as.factor(batch_vector)

# 计算离群值 (距离中心的欧氏距离)
center <- colMeans(pca_df[, c("PC1", "PC2")])
distances <- sqrt(rowSums((pca_df[, c("PC1", "PC2")] - center)^2))
pca_df$Distance <- distances

# 标记离群样本 (距离最大的前5%)
outlier_threshold <- quantile(distances, 0.95)
pca_df$Outlier <- ifelse(distances > outlier_threshold, "Yes", "No")


# 创建带标签的PCA图，无batch图注和网格
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Dataset)) +
  geom_point(size = 3, alpha = 0.8, aes(shape = Batch)) +
  geom_label_repel(
    data = subset(pca_df, Outlier == "Yes"),
    aes(label = Sample),
    box.padding = 0.5,
    segment.color = "grey50",
    min.segment.length = 0,
    max.overlaps = 20,
    size = 2.5,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = 1:length(unique(batch_vector)), guide = "none") +  # 隐藏Batch图注
  labs(
    title = "PCA of Batch-Corrected Expression Data",
    subtitle = paste("Labeled outliers (top 5% distance from center, threshold =", 
                     round(outlier_threshold, 2)),
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank()   # 移除次要网格线
  ) +
  guides(color = guide_legend(nrow = 2))  # 只保留Dataset图注

pca_plot














# 保存PCA图为PDF
ggsave("PCA_5sample.pdf", pca_plot, width = 13, height = 11)


# 8. 识别并输出离群样本
outliers <- pca_df %>% 
  filter(Outlier == "Yes") %>% 
  arrange(desc(Distance))

cat("Identified outliers (potential samples to remove):\n")
print(outliers[, c("Sample", "Dataset", "Batch", "Distance")])

# 9. 剔除离群样本 (可选)
clean_expression <- final_expression[!rownames(final_expression) %in% outliers$Sample, ]

dim(clean_expression)

# 保存去批次且去除离群值后的矩阵
write.csv(clean_expression, "batch_corrected_expression_clean.csv")

#####################################################################################
#####################################################################################
#####################################################################################



# 读取去除离群值后的数据
clean_expression <- read.csv("batch_corrected_expression_clean.csv", row.names = 1)

# 进行PCA分析
pca_result <- prcomp(clean_expression, scale. = TRUE)

# 创建PCA结果数据框
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Sample <- rownames(pca_df)
pca_df$Dataset <- sapply(strsplit(pca_df$Sample, "_"), `[`, 1)

# 重新创建批次映射
batch_mapping <- c(
  "GSE7670" = 1,
  "GSE18842" = 2,
  "GSE19804" = 3,
  "GSE31547" = 4,
  "GSE33532" = 5
)
pca_df$Batch <- as.factor(batch_mapping[pca_df$Dataset])

# 计算离群值 (距离中心的欧氏距离) - 注意这里仅用于绘图逻辑，实际已无离群值
center <- colMeans(pca_df[, c("PC1", "PC2")])
distances <- sqrt(rowSums((pca_df[, c("PC1", "PC2")] - center)^2))
pca_df$Distance <- distances
outlier_threshold <- quantile(distances, 0.95)  # 保留阈值计算用于标题
pca_df$Outlier <- "No"  # 所有样本标记为非离群值

# 使用原始绘图代码绘制PCA图
library(ggplot2)
library(ggrepel)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Dataset)) +
  geom_point(size = 3, alpha = 0.8, aes(shape = Batch)) +
  geom_label_repel(
    data = subset(pca_df, Outlier == "Yes"),  # 这里不会有任何样本
    aes(label = Sample),
    box.padding = 0.5,
    segment.color = "grey50",
    min.segment.length = 0,
    max.overlaps = 20,
    size = 2.5,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = 1:5, guide = "none") +  # 隐藏Batch图注
  labs(
    title = "PCA of Batch-Corrected Expression Data (Outliers Removed)",
    subtitle = paste("Outlier removal threshold =", round(outlier_threshold, 2)),
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank()   # 移除次要网格线
  ) +
  guides(color = guide_legend(nrow = 2))  # 只保留Dataset图注

# 显示图形
pca_plot

# 保存为PDF
ggsave("PCA_clean_5sample.pdf", pca_plot, width = 13, height = 11)

