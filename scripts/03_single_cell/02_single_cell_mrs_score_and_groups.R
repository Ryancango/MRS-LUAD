# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/单细胞风险评分调整spock2和adh1b.R
# Role: Calculate single-cell MRS scores using the final coefficient set and assign median-based risk groups.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}

library(Seurat)

# Final coefficient set used in the manuscript and downstream single-cell scoring.
if (!exists("mrs_coefficients")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}
genes <- names(mrs_coefficients)
average_coef <- mrs_coefficients[genes]

#需要预先加载数据
# 假设已加载必要的包，如library(Seurat)
# Example: load the processed GSE131907 lung Seurat object as seurat_tlung_riskgroup before running.
for (obj_name in c("seurat_tlung_riskgroup")) {
  obj <- get(obj_name)
  expr_matrix <- GetAssayData(obj, assay = "RNA", layer = "data") # 规范化表达
  missing_genes <- setdiff(genes, rownames(expr_matrix))
  if (length(missing_genes) > 0) {
    cat(paste("Missing genes in", obj_name, ":", paste(missing_genes, collapse = ", "), "\n"))
    # 填充0
    for (g in missing_genes) {
      expr_matrix <- rbind(expr_matrix, rep(0, ncol(expr_matrix)))
      rownames(expr_matrix)[nrow(expr_matrix)] <- g
    }
  }
  expr_genes <- expr_matrix[genes, ]
  obj$risk_score <- colSums(expr_genes * average_coef[genes]) # 计算risk_score
  
  
  # 使用50%分位数作为阈值
  threshold_risk <- quantile(obj$risk_score, 0.5, na.rm = TRUE)
  obj$risk_group <- ifelse(obj$risk_score > threshold_risk, "High Risk", "Low Risk")
  assign(obj_name, obj) # 更新对象
  
  # 使用50%分位数作为阈值
  #median_risk <- median(obj$risk_score, na.rm = TRUE)
  #obj$risk_group <- ifelse(obj$risk_score > median_risk, "High Risk", "Low Risk")
  #assign(obj_name, obj) # 更新对象
  
  # 可视化tSNE，按risk_group着色，并指定颜色
  # 定义颜色：Low Risk 用 #8CDEBE，High Risk 用 tomato
  custom_cols <- c("Low Risk" = "#8CDEBE", "High Risk" = "tomato")
  
  p <- DimPlot(obj, reduction = "tsne", group.by = "risk_group", label = TRUE, pt.size = 0.3, cols = custom_cols)
  print(p)
  #ggsave(paste0(obj_name, "_tsne_risk_group.png"), plot = p, width = 8, height = 6) # 保存图片
}



######risk_group在每个Cell_type中的占比
######risk_group在每个Cell_type中的占比
######risk_group在每个Cell_type中的占比

risk_obj <- seurat_tlung_riskgroup
cell_type_risk <- table(risk_obj$Cell_type, risk_obj$risk_group)
prop_table <- prop.table(cell_type_risk, margin = 1) * 100  # 按Cell_type计算百分比
print(prop_table)
library(ggplot2)
prop_df <- as.data.frame(prop_table)
colnames(prop_df) <- c("Cell_type", "risk_group", "Percentage")
p <- ggplot(prop_df, aes(x = Cell_type, y = Percentage, fill = risk_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("High Risk" = "tomato", "Low Risk" = "#8CDEBE")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "Percentage (%)", fill = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)







