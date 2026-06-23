# Curated release script for: MRS-LUAD
# Source script: intersection_mechine/mechine/outside_data/mime1**.R
# Role: Train and compare 117 machine-learning survival models using TCGA-LUAD, GSE31210 and GSE72094 splits.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}

library(Mime1)
# 读取候选基因CSV
genelist <- read.csv("significant_genes_original_alias.csv")$gene

load(c(TCGALUAD.rdata,GSE72094.RData,GSE31210.RData))

# 计算共同候选基因（确保所有数据集都有）
common_genes <- Reduce(intersect, list(intersect(genelist, colnames(tcgaluad)),
                                       intersect(genelist, colnames(GSE31210)),
                                       intersect(genelist, colnames(GSE72094))))

# 计算共同候选基因（确保所有数据集都有）
common_genes <- Reduce(intersect, list(intersect(genelist, colnames(tcgaluad)),
                                       intersect(genelist, colnames(GSE72094)),
                                       intersect(genelist, colnames(GSE31210))))

# 修改filter_data函数，使用common_genes，并过滤OS.time > 0
filter_data <- function(df, genes) {
  survival_cols <- c("ID", "OS.time", "OS")
  gene_cols <- intersect(genes, colnames(df))
  df <- df[df$OS.time > 0, ]  # 移除OS.time <= 0的样本
  df[, c(survival_cols, gene_cols)]
}

tcgaluad_filtered <- filter_data(tcgaluad, common_genes)
GSE31210_filtered <- filter_data(GSE31210, common_genes)
GSE72094_filtered <- filter_data(GSE72094, common_genes)

# 分割每个数据集（50/50）
set.seed(5201314)
train_idx_tcga <- sample(nrow(tcgaluad_filtered), size = round(0.5 * nrow(tcgaluad_filtered)))
tcga_train <- tcgaluad_filtered[train_idx_tcga, ]
tcga_vali <- tcgaluad_filtered[-train_idx_tcga, ]

train_idx_gse31210 <- sample(nrow(GSE31210_filtered), size = round(0.5 * nrow(GSE31210_filtered)))
gse31210_train <- GSE31210_filtered[train_idx_gse31210, ]
gse31210_vali <- GSE31210_filtered[-train_idx_gse31210, ]

train_idx_gse72094 <- sample(nrow(GSE72094_filtered), size = round(0.5 * nrow(GSE72094_filtered)))
gse72094_train <- GSE72094_filtered[train_idx_gse72094, ]
gse72094_vali <- GSE72094_filtered[-train_idx_gse72094, ]


# 合并验证部分
vali_data <- rbind(tcga_vali, gse31210_vali, gse72094_vali)
train_data <- rbind(tcga_train, gse31210_train, gse72094_train)

# 创建数据集列表（包括合并train和各vali）
list_train_vali_Data <- list(
  #train_data = train_data,
  TCGA_LUAD_Vali = tcga_vali,
  GSE31210_Vali = gse31210_vali,
  GSE72094_Vali = gse72094_vali
  #tcgalusc_Vali = tcgalusc_vali
)

# 构建预后模型（其余不变）
res55 <- ML.Dev.Prog.Sig(
  train_data = train_data,
  list_train_vali_Data = list_train_vali_Data,
  unicox.filter.for.candi = TRUE,
  unicox_p_cutoff = 0.05,
  candidate_genes = common_genes,  # 更新为共同基因
  mode = 'all',  # 'all'为所有组合，'single'或'double'可选
  nodesize = 5,
  seed = 1314  # 固定种子以复现
)

# 可视化C-index分布（调整validate_set和order）
cindex_dis_all(res55, validate_set = c("TCGA_LUAD_Vali",
                                       "GSE31210_Vali", 
                                       "GSE72094_Vali",
                                       "train_data"
                                       ), 
               order = names(list_train_vali_Data), width = 0.35)


#筛选核心基因
res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = vali_data,
                                              candidate_genes = genelist,
                                              mode = "all",nodesize =5,seed = 1314 )

p3 <- core_feature_rank(res.feature.all, top=20)
ggsave(filename = "gene55.pdf", plot = p3, width = 10, height = 8, dpi = 300)






# 计算1年AUC并可视化（调整dataset, validate_set, order；类似可改3年、5年）
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res, train_data = train_data,
                             inputmatrix.list = list_train_vali_Data, mode = 'all', AUC_time = 1, auc_cal_method = "KM")
auc_dis_all(all.auc.1y, dataset = names(list_train_vali_Data), 
            validate_set = c("TCGA_LUAD_Vali", "GSE31210_Vali", "tcgalusc_Vali"),
            order = names(list_train_vali_Data), width = 0.35, year = 1)


# 示例：绘制最佳模型的生存曲线（假设最佳模型为"StepCox[forward] + plsRcox"，运行后根据res调整；调整dataset names）
survplot <- vector("list", length(vali_data))
names(survplot) <- names(vali_data)
for (i in seq_along(vali_data)) {
  survplot[[i]] <- rs_sur(res55, model_name = "CoxBoost + survival−SVM", dataset = names(vali_data)[i],
                          median.line = "hv", cutoff = 0.5, conf.int = TRUE, xlab = "Day", pval.coord = c(1000, 0.9))
}
aplot::plot_list(gglist = survplot, ncol = 2)


