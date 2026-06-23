# Curated release script for: MRS-LUAD
# Source script: scRNA/131907/Pseudotime_Quantification/monocle2.R
# Role: Infer epithelial-cell pseudotime trajectories with Monocle2.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}

# 重新分析：关注轨迹状态、进化时序和组织来源
library(Seurat)
library(monocle)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)
library(tidyverse)


# 设置SCI级别的绘图主题 - 更大的字体和元素，无网格线
sci_theme <- theme_bw(base_size = 17) +  # 基础字体放大3号
  theme(
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.text = element_text(color = "black", size = 15),  # 放大3号
    axis.title = element_text(color = "black", size = 16, face = "bold"),  # 放大3号
    legend.title = element_text(size = 15, face = "bold"),  # 放大3号
    legend.text = element_text(size = 14),  # 放大3号
    legend.key.size = unit(0.6, "cm"),
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),  # 放大3号
    plot.subtitle = element_text(size = 15, hjust = 0.5),  # 放大3号
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.7),
    strip.text = element_text(size = 14)  # 放大3号
  )

# 1. 加载数据
cat("=== 加载数据 ===\n")
seurat_obj <- readRDS("seurat_tlungepi_risk.rds")
#seurat_obj <- readRDS("seurat_lung_risk.rds")

# 检查数据基本信息
cat("数据维度:", dim(seurat_obj), "\n")
cat("患者来源分布:\n")
print(table(seurat_obj$Sample_Origin))
cat("风险组分布:\n")
print(table(seurat_obj$risk_group))

# 2. 使用所有细胞进行拟时序分析
cat("\n=== 对所有细胞进行拟时序分析 ===\n")

# 提取表达矩阵
expression_matrix <- as.matrix(seurat_obj@assays$RNA@layers$counts)
rownames(expression_matrix) <- rownames(seurat_obj)
colnames(expression_matrix) <- colnames(seurat_obj)

# 创建CellDataSet对象
pd <- new("AnnotatedDataFrame", data = seurat_obj@meta.data)
fData <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(fData) <- rownames(expression_matrix)
fd <- new("AnnotatedDataFrame", data = fData)

cds <- newCellDataSet(expression_matrix,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily = negbinomial.size())

# Monocle2预处理
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 选择高变基因 - 使用更严格的阈值以获得更好的轨迹
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, 
                         mean_expression >= 0.5 &  # 提高表达量阈值
                           dispersion_empirical >= 1.5 * dispersion_fit)$gene_id  # 提高离散度阈值

cds <- setOrderingFilter(cds, ordering_genes)
cat("用于轨迹分析的高变基因数量:", length(ordering_genes), "\n")

# 降维 - 使用DDRTree方法
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree', norm_method = 'log')

# 3. 智能选择根状态
cat("\n=== 智能选择根状态 ===\n")

# 首先运行一次不指定根状态的orderCells来设置状态
cds <- orderCells(cds)

# 分析每个状态的特征
state_analysis <- data.frame(
  State = cds$State,
  Sample_Origin = seurat_obj$Sample_Origin,
  Risk_Score = seurat_obj$risk_score
) %>%
  group_by(State) %>%
  summarise(
    n_cells = n(),
    tLung_prop = sum(Sample_Origin == "tLung") / n() * 100,
    nLung_prop = sum(Sample_Origin == "nLung") / n() * 100,
    mean_risk = mean(Risk_Score)
  )

cat("各状态特征分析:\n")
print(state_analysis)

# 基于生物学知识选择根状态：
# 通常选择正常组织细胞比例高的状态作为起点
if(max(state_analysis$nLung_prop) > 50) {
  # 如果有明显正常细胞富集的状态，选择它作为根
  root_state <- state_analysis$State[which.max(state_analysis$nLung_prop)]
  cat("选择正常细胞富集的状态作为根状态: State", root_state, "\n")
} else {
  # 否则选择细胞数量最多的状态
  root_state <- state_analysis$State[which.max(state_analysis$n_cells)]
  cat("选择细胞数量最多的状态作为根状态: State", root_state, "\n")
}

# 使用选择的根状态重新排序细胞
cds <- orderCells(cds, root_state = root_state)

# 4. 提取轨迹数据并分析
cat("\n=== 提取轨迹数据 ===\n")
ddrtree_coords <- reducedDimS(cds)
pseudotime_values <- cds$Pseudotime

# 创建轨迹数据框
ddrtree_df <- data.frame(
  Component1 = ddrtree_coords[1, ],
  Component2 = ddrtree_coords[2, ],
  Pseudotime = pseudotime_values,
  Risk_Score = seurat_obj$risk_score,
  Risk_Group = seurat_obj$risk_group,
  Sample_Origin = seurat_obj$Sample_Origin,
  State = as.factor(cds$State),
  row.names = colnames(cds)
)

# 5. 计算各状态的细胞分布
cat("\n=== 计算各状态的细胞分布 ===\n")

state_summary <- ddrtree_df %>%
  group_by(State) %>%
  summarise(
    n_cells = n(),
    tLung_Proportion = sum(Sample_Origin == "tLung") / n() * 100,
    HighRisk_Proportion = sum(Risk_Group == "High Risk") / n() * 100,
    Mean_Risk_Score = mean(Risk_Score),
    Mean_Pseudotime = mean(Pseudotime)
  ) %>%
  arrange(Mean_Pseudotime)  # 按平均拟时序排序

cat("按拟时序排序的各状态分布:\n")
print(state_summary)

# 6. 创建高质量的可视化 - 更大的图片
cat("\n=== 创建高质量可视化 ===\n")

# 设置更大的图片尺寸
single_width <- 10  # 增加宽度
single_height <- 8  # 增加高度
combined_width <- 18
combined_height <- 14

# 颜色方案
tumor_color <- "#D62728"  # 红色代表肿瘤
normal_color <- "#1F77B4"  # 蓝色代表正常
high_risk_color <- "#FF7F0E"  # 橙色代表高风险
low_risk_color <- "#2CA02C"   # 绿色代表低风险

# 图1: 轨迹图 - 按患者来源着色
p1 <- ggplot(ddrtree_df, aes(x = Component1, y = Component2, color = Sample_Origin)) +
  geom_point(size = 1.8, alpha = 0.7) +  # 增大点的大小
  scale_color_manual(values = c("tLung" = tumor_color, "nLung" = normal_color),
                     name = "Tissue Origin",
                     labels = c("tLung" = "Tumor", "nLung" = "Normal")) +
  labs(title = "Cell Trajectory by Tissue Origin",
       x = "Component 1", y = "Component 2") +
  sci_theme

# 图2: 轨迹图 - 按状态着色，添加状态标签
state_centers <- ddrtree_df %>%
  group_by(State) %>%
  summarise(
    center_x = median(Component1),
    center_y = median(Component2)
  )

p2 <- ggplot(ddrtree_df, aes(x = Component1, y = Component2, color = State)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_label(data = state_centers, 
             aes(x = center_x, y = center_y, label = paste("State", State)),
             size = 6, fontface = "bold", alpha = 0.8) +  # 增大标签字体
  scale_color_manual(values = c("1" = "#4DAF4A", "2" = "#FF7F00", "3" = "#E41A1C")) +
  labs(title = "Cell Trajectory by State",
       x = "Component 1", y = "Component 2") +
  sci_theme

# 图3: 轨迹图 - 按拟时序着色，添加进化方向箭头
p3 <- ggplot(ddrtree_df, aes(x = Component1, y = Component2, color = Pseudotime)) +
  geom_point(size = 1.8, alpha = 0.7) +
  scale_color_viridis_c(name = "Pseudotime", direction = -1) +
  # 添加进化方向箭头
  annotate("segment", 
           x = min(ddrtree_df$Component1), 
           xend = min(ddrtree_df$Component1) + diff(range(ddrtree_df$Component1)) * 0.15,
           y = min(ddrtree_df$Component2),
           yend = min(ddrtree_df$Component2),
           arrow = arrow(type = "closed", length = unit(0.3, "cm")),
           color = "black", size = 1.5) +
  annotate("text", 
           x = min(ddrtree_df$Component1) + diff(range(ddrtree_df$Component1)) * 0.08,
           y = min(ddrtree_df$Component2) - diff(range(ddrtree_df$Component2)) * 0.05,
           label = "Developmental\ndirection", 
           size = 5, fontface = "bold") +  # 增大箭头文字
  labs(title = "Developmental Trajectory by Pseudotime",
       x = "Component 1", y = "Component 2") +
  sci_theme

# 新图: 轨迹图 - 按风险评分着色，高风险点在最上层
cat("\n=== 创建风险评分轨迹图 ===\n")

# 按风险评分排序，确保高风险点最后绘制（在最上层）
ddrtree_df_sorted <- ddrtree_df[order(ddrtree_df$Risk_Score), ]

p_risk <- ggplot(ddrtree_df_sorted, aes(x = Component1, y = Component2, color = Risk_Score)) +
  geom_point(size = 1.8, alpha = 0.7) +
  scale_color_gradient2(
    low = low_risk_color,
    mid = "white",
    high = high_risk_color,
    midpoint = median(ddrtree_df$Risk_Score),
    name = "Risk Score"
  ) +
  labs(title = "Cell Trajectory by Risk Score",
       subtitle = "High-risk cells displayed on top",
       x = "Component 1", y = "Component 2") +
  sci_theme

# 图4: 状态比例对比图
plot_data <- state_summary %>%
  select(State, tLung_Proportion, HighRisk_Proportion) %>%
  pivot_longer(cols = c(tLung_Proportion, HighRisk_Proportion), 
               names_to = "Type", values_to = "Proportion")

p4 <- ggplot(plot_data, aes(x = State, y = Proportion, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Proportion)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 5.5, fontface = "bold") +  # 增大文字
  scale_fill_manual(values = c("tLung_Proportion" = tumor_color, 
                               "HighRisk_Proportion" = high_risk_color),
                    labels = c("tLung_Proportion" = "Tumor Cells", 
                               "HighRisk_Proportion" = "High-Risk Cells")) +
  labs(title = "Cell Type Proportions by Trajectory State",
       x = "Trajectory State", y = "Proportion (%)") +
  ylim(0, max(plot_data$Proportion) * 1.15) +
  sci_theme

# 图5: 拟时序与风险评分的关系，按组织来源分面
p5 <- ggplot(ddrtree_df, aes(x = Pseudotime, y = Risk_Score)) +
  geom_point(aes(color = Sample_Origin), alpha = 0.6, size = 1.5) +  # 增大点
  geom_smooth(aes(color = Sample_Origin, fill = Sample_Origin), 
              method = "loess", se = TRUE, alpha = 0.2, linewidth = 1.2) +
  facet_wrap(~Sample_Origin, ncol = 2) +
  scale_color_manual(values = c("tLung" = tumor_color, "nLung" = normal_color),
                     name = "Tissue Origin",
                     labels = c("tLung" = "Tumor", "nLung" = "Normal")) +
  scale_fill_manual(values = c("tLung" = tumor_color, "nLung" = normal_color),
                    name = "Tissue Origin",
                    labels = c("tLung" = "Tumor", "nLung" = "Normal")) +
  labs(title = "Risk Score Dynamics along Developmental Trajectory",
       x = "Pseudotime", y = "Risk Score") +
  sci_theme

# 图6: 拟时序在各状态的分布
p6 <- ggplot(ddrtree_df, aes(x = State, y = Pseudotime, fill = State)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red") +  # 增大点
  labs(title = "Pseudotime Distribution by State",
       x = "Trajectory State", y = "Pseudotime",
       subtitle = "Red diamond indicates mean pseudotime") +
  sci_theme +
  theme(legend.position = "none")

# 7. 相关系数分析
cat("\n=== 进行相关系数分析 ===\n")

# 创建二元变量用于相关性分析
ddrtree_df$tLung_binary <- ifelse(ddrtree_df$Sample_Origin == "tLung", 1, 0)
ddrtree_df$HighRisk_binary <- ifelse(ddrtree_df$Risk_Group == "High Risk", 1, 0)

# 计算tLung与高风险的相关性
cor_test_tLung_highrisk <- cor.test(ddrtree_df$tLung_binary, ddrtree_df$HighRisk_binary)

# 计算与拟时序的相关性
cor_test_tLung_pseudotime <- cor.test(ddrtree_df$tLung_binary, ddrtree_df$Pseudotime)
cor_test_highrisk_pseudotime <- cor.test(ddrtree_df$HighRisk_binary, ddrtree_df$Pseudotime)

# 创建相关性结果数据框
correlation_results <- data.frame(
  Comparison = c("tLung vs HighRisk", "tLung vs Pseudotime", "HighRisk vs Pseudotime"),
  Correlation_Coefficient = c(
    cor_test_tLung_highrisk$estimate,
    cor_test_tLung_pseudotime$estimate,
    cor_test_highrisk_pseudotime$estimate
  ),
  P_value = c(
    cor_test_tLung_highrisk$p.value,
    cor_test_tLung_pseudotime$p.value,
    cor_test_highrisk_pseudotime$p.value
  ),
  Method = c(
    as.character(cor_test_tLung_highrisk$method),
    as.character(cor_test_tLung_pseudotime$method),
    as.character(cor_test_highrisk_pseudotime$method)
  )
)

# 添加显著性标记
correlation_results$Significance <- ifelse(correlation_results$P_value < 0.001, "***",
                                           ifelse(correlation_results$P_value < 0.01, "**",
                                                  ifelse(correlation_results$P_value < 0.05, "*", "NS")))

cat("相关性分析结果:\n")
print(correlation_results)

# 创建相关性可视化图
p_correlation <- ggplot(correlation_results, aes(x = Comparison, y = Correlation_Coefficient)) +
  geom_col(aes(fill = Correlation_Coefficient), alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("r = %.3f\np = %.2e\n%s", 
                                Correlation_Coefficient, P_value, Significance)),
            vjust = -0.3, size = 5.5, fontface = "bold", lineheight = 0.8) +  # 增大文字
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Correlation") +
  labs(title = "Correlation Analysis: tLung, HighRisk and Pseudotime",
       x = "Variable Pairs", y = "Correlation Coefficient") +
  ylim(min(correlation_results$Correlation_Coefficient) * 1.2, 
       max(correlation_results$Correlation_Coefficient) * 1.2) +
  sci_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 创建散点图矩阵展示关系
p_scatter_matrix <- ggplot(ddrtree_df, aes(x = tLung_binary, y = HighRisk_binary)) +
  geom_jitter(aes(color = Pseudotime), alpha = 0.6, size = 2, width = 0.1, height = 0.1) +  # 增大点
  scale_color_viridis_c(name = "Pseudotime") +
  labs(title = "Relationship between tLung and HighRisk Status",
       subtitle = sprintf("Correlation: r = %.3f, p = %.2e", 
                          cor_test_tLung_highrisk$estimate, 
                          cor_test_tLung_highrisk$p.value),
       x = "tLung Status (0 = nLung, 1 = tLung)", 
       y = "HighRisk Status (0 = Low Risk, 1 = High Risk)") +
  sci_theme

# 8. 保存高质量PDF图片
cat("\n=== 保存高质量PDF图片 ===\n")

output_dir <- "test"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 保存单个图片
pdf(file.path(output_dir, "Figure1_Trajectory_TissueOrigin.pdf"), 
    width = single_width, height = single_height)
print(p1)
dev.off()
cat("保存: Figure1_Trajectory_TissueOrigin.pdf\n")

pdf(file.path(output_dir, "Figure2_Trajectory_State.pdf"), 
    width = single_width, height = single_height)
print(p2)
dev.off()
cat("保存: Figure2_Trajectory_State.pdf\n")

pdf(file.path(output_dir, "Figure3_Trajectory_Pseudotime.pdf"), 
    width = single_width, height = single_height)
print(p3)
dev.off()
cat("保存: Figure3_Trajectory_Pseudotime.pdf\n")

pdf(file.path(output_dir, "Figure4_Trajectory_RiskScore.pdf"), 
    width = single_width, height = single_height)
print(p_risk)
dev.off()
cat("保存: Figure4_Trajectory_RiskScore.pdf\n")

pdf(file.path(output_dir, "Figure5_State_Proportions.pdf"), 
    width = single_width, height = single_height)
print(p4)
dev.off()
cat("保存: Figure5_State_Proportions.pdf\n")

pdf(file.path(output_dir, "Figure6_RiskScore_Dynamics.pdf"), 
    width = single_width, height = single_height)
print(p5)
dev.off()
cat("保存: Figure6_RiskScore_Dynamics.pdf\n")

pdf(file.path(output_dir, "Figure7_Pseudotime_Distribution.pdf"), 
    width = single_width, height = single_height)
print(p6)
dev.off()
cat("保存: Figure7_Pseudotime_Distribution.pdf\n")

pdf(file.path(output_dir, "Figure8_Correlation_Analysis.pdf"), 
    width = single_width, height = single_height)
print(p_correlation)
dev.off()
cat("保存: Figure8_Correlation_Analysis.pdf\n")

pdf(file.path(output_dir, "Figure9_Scatter_Matrix.pdf"), 
    width = single_width, height = single_height)
print(p_scatter_matrix)
dev.off()
cat("保存: Figure9_Scatter_Matrix.pdf\n")

# 组合图1: 主要轨迹分析
combined_main <- (p1 + p2) / (p3 + p_risk) +
  plot_annotation(tag_levels = 'A', 
                  title = 'Comprehensive Analysis of Lung Epithelial Cell Trajectory',
                  theme = theme(plot.title = element_text(size = 21, face = "bold", hjust = 0.5)))  # 增大标题

pdf(file.path(output_dir, "Figure10_Combined_Main.pdf"), 
    width = combined_width, height = combined_height)
print(combined_main)
dev.off()
cat("保存: Figure10_Combined_Main.pdf\n")

# 组合图2: 动态分析
combined_dynamics <- (p5 + p6) +
  plot_annotation(tag_levels = 'A',
                  title = 'Temporal Dynamics and Risk Assessment',
                  theme = theme(plot.title = element_text(size = 21, face = "bold", hjust = 0.5)))

pdf(file.path(output_dir, "Figure11_Combined_Dynamics.pdf"), 
    width = combined_width, height = combined_height/1.5)
print(combined_dynamics)
dev.off()
cat("保存: Figure11_Combined_Dynamics.pdf\n")

# 组合图3: 相关性分析
combined_correlation <- (p_correlation + p_scatter_matrix) +
  plot_annotation(tag_levels = 'A',
                  title = 'Correlation Analysis and Relationship Mapping',
                  theme = theme(plot.title = element_text(size = 21, face = "bold", hjust = 0.5)))

pdf(file.path(output_dir, "Figure12_Combined_Correlation.pdf"), 
    width = combined_width, height = combined_height/1.8)
print(combined_correlation)
dev.off()
cat("保存: Figure12_Combined_Correlation.pdf\n")

# 9. 保存分析结果
write.csv(state_summary, file.path(output_dir, "State_Analysis_Summary.csv"), row.names = FALSE)
write.csv(correlation_results, file.path(output_dir, "Correlation_Analysis_Results.csv"), row.names = FALSE)
saveRDS(cds, file.path(output_dir, "monocle_cds_optimized.rds"))
saveRDS(ddrtree_df, file.path(output_dir, "trajectory_data_optimized.rds"))

# 10. 生成分析报告
sink(file.path(output_dir, "Analysis_Report.txt"))
cat("=== REVISED PSEUDOTIME ANALYSIS REPORT ===\n\n")
cat("Analysis Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("Total cells analyzed:", ncol(seurat_obj), "\n")
cat("High-variable genes used:", length(ordering_genes), "\n")
cat("Root state selected:", root_state, "\n\n")

cat("=== CORRELATION ANALYSIS RESULTS ===\n")
print(correlation_results)
cat("\n")

cat("=== KEY CORRELATION INSIGHTS ===\n")
cat(sprintf("1. tLung vs HighRisk: r = %.3f, p = %.2e\n", 
            correlation_results$Correlation_Coefficient[1],
            correlation_results$P_value[1]))
cat(sprintf("2. tLung vs Pseudotime: r = %.3f, p = %.2e\n", 
            correlation_results$Correlation_Coefficient[2],
            correlation_results$P_value[2]))
cat(sprintf("3. HighRisk vs Pseudotime: r = %.3f, p = %.2e\n", 
            correlation_results$Correlation_Coefficient[3],
            correlation_results$P_value[3]))
cat("\n")

cat("=== BIOLOGICAL INTERPRETATION ===\n")
if(correlation_results$P_value[1] < 0.05) {
  if(correlation_results$Correlation_Coefficient[1] > 0) {
    cat("- tLung and HighRisk show significant positive correlation, suggesting tumor cells are enriched in high-risk population.\n")
  } else {
    cat("- tLung and HighRisk show significant negative correlation.\n")
  }
} else {
  cat("- No significant correlation between tLung and HighRisk status.\n")
}

if(correlation_results$P_value[2] < 0.05) {
  if(correlation_results$Correlation_Coefficient[2] > 0) {
    cat("- tLung cells are significantly enriched along pseudotime progression.\n")
  } else {
    cat("- tLung cells are significantly depleted along pseudotime progression.\n")
  }
}

if(correlation_results$P_value[3] < 0.05) {
  if(correlation_results$Correlation_Coefficient[3] > 0) {
    cat("- HighRisk cells are significantly enriched along pseudotime progression.\n")
  } else {
    cat("- HighRisk cells are significantly depleted along pseudotime progression.\n")
  }
}
cat("\n")

cat("=== TRAJECTORY STATE CHARACTERISTICS ===\n")
print(state_summary)
cat("\n")

cat("=== KEY BIOLOGICAL INSIGHTS ===\n")

# 分析状态间的进化关系
early_state <- state_summary$State[which.min(state_summary$Mean_Pseudotime)]
late_state <- state_summary$State[which.max(state_summary$Mean_Pseudotime)]

cat("Developmental progression: State", early_state, "-> State", late_state, "\n")
cat(sprintf("State %d (early): %.1f%% tumor cells, %.1f%% high-risk cells\n", 
            early_state,
            state_summary$tLung_Proportion[state_summary$State == early_state],
            state_summary$HighRisk_Proportion[state_summary$State == early_state]))
cat(sprintf("State %d (late): %.1f%% tumor cells, %.1f%% high-risk cells\n", 
            late_state,
            state_summary$tLung_Proportion[state_summary$State == late_state],
            state_summary$HighRisk_Proportion[state_summary$State == late_state]))

# 计算变化趋势
tumor_change <- state_summary$tLung_Proportion[state_summary$State == late_state] - 
  state_summary$tLung_Proportion[state_summary$State == early_state]
risk_change <- state_summary$HighRisk_Proportion[state_summary$State == late_state] - 
  state_summary$HighRisk_Proportion[state_summary$State == early_state]

cat(sprintf("Tumor cell enrichment along trajectory: %.1f%% increase\n", tumor_change))
cat(sprintf("High-risk cell enrichment along trajectory: %.1f%% increase\n", risk_change))

sink()

cat("\n=== 分析完成 ===\n")
cat("输出目录:", output_dir, "\n")
cat("生成的图片:\n")
cat("- 9个单图 (Figure1-9_*.pdf)\n")
cat("- 3个组合图 (Figure10-12_*.pdf)\n")
cat("\n关键改进:\n")
cat("- 新增风险评分轨迹图，高风险点显示在最上层\n")
cat("- 所有文字放大3号，适合SCI发表\n")
cat("- 去掉背景网格线，提升图片清晰度\n")
cat("- 增大图片尺寸，提高可读性\n")
cat("- 添加了完整的相关性分析\n")
cat("- 重新整理和优化代码结构\n")