# Curated release script for: MRS-LUAD
# Source script: intersection_mechine/mechine/outside_data/xena/pancancerKMROC.R
# Role: Assess 10-gene panel prognostic associations across TCGA cancer types and endpoints.
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
# 0. 加载必要的 R 包
# ==============================================================================
library(survival)
library(survminer)
library(ggplot2)
library(survivalROC)
library(dplyr)
library(grid)

# ==============================================================================
# 1. 基础设置 (基因、生存类型、输入输出)
# ==============================================================================

# A. 设置工作目录 (可选，建议设置到你的数据所在目录)
# setwd("/你的/工作/目录/")

# B. 定义基因列表 (你的 10 个标志物)
genes <- c("RRAS", "GPM6B", "FMO3", "CHRDL1", "ADRB2", 
           "LPL", "ANGPT1", "SPOCK2", "P2RY14", "ADH1B")

# C. 定义 4 种生存数据类型
surv_types <- list(
  "OS"  = c("OS.time", "OS"),
  "DSS" = c("DSS.time", "DSS"),
  "DFI" = c("DFI.time", "DFI"),
  "PFI" = c("PFI.time", "PFI")
)

# D. 定义输出目录名称
output_dir <- "01_SCI_Plots_Gallery"

# 如果目录不存在，创建它
if (!dir.exists(output_dir)) dir.create(output_dir)

# ==============================================================================
# 2. 读取数据
# ==============================================================================
# 请修改下面的路径为你 RDS 文件的实际路径
rds_path <- file.path(mrs_project_root(), "intersection_mechine/mechine/outside_data/xena/TCGA_PANCAN_Merged_Final.rds")

if (file.exists(rds_path)) {
  cat("正在读取数据文件:", rds_path, "...\n")
  TCGA_PANCAN_Merged_Final <- readRDS(rds_path)
  cat("数据读取成功！包含样本数:", nrow(TCGA_PANCAN_Merged_Final), "\n")
} else {
  stop("错误：找不到 RDS 文件，请检查路径是否正确！")
}

# ==============================================================================
# 3. 定义核心绘图函数 (高度封装)
# ==============================================================================

analyze_and_plot <- function(data, cancer_name, type_name, time_col, status_col, gene_list, save_path) {
  
  # --- 步骤 1: 数据清洗与检查 ---
  # 提取相关列
  req_cols <- c(time_col, status_col, gene_list)
  # 检查列是否存在
  if (!all(req_cols %in% colnames(data))) return(NULL)
  
  # 去除 NA
  sub_data <- na.omit(data[, req_cols])
  
  # 检查样本量 (至少 20 个样本，且至少 5 个死亡/事件发生，否则做 Cox 回归不可靠)
  if (nrow(sub_data) < 20 || sum(sub_data[[status_col]] == 1) < 5) return(NULL)
  
  # 检查基因方差 (剔除在当前癌症中不表达或表达量无变化的基因)
  valid_genes <- gene_list[sapply(sub_data[, gene_list], var) > 0]
  if (length(valid_genes) == 0) return(NULL)
  
  # --- 步骤 2: 构建多因素 Cox 模型并计算 Risk Score ---
  surv_obj <- paste0("Surv(", time_col, ", ", status_col, ")")
  form_cox <- as.formula(paste(surv_obj, "~", paste(valid_genes, collapse = " + ")))
  
  # 使用 tryCatch 防止模型不收敛导致报错
  cox_model <- tryCatch({
    coxph(form_cox, data = sub_data, control = coxph.control(iter.max = 100))
  }, error = function(e) NULL, warning = function(w) NULL)
  
  if (is.null(cox_model) || any(is.na(coef(cox_model)))) return(NULL)
  
  # 计算 Risk Score
  coefficients <- coef(cox_model)
  # 矩阵乘法计算得分
  risk_scores <- as.matrix(sub_data[, valid_genes]) %*% coefficients
  sub_data$risk_score <- as.vector(risk_scores)
  
  # 根据中位数分组
  median_risk <- median(sub_data$risk_score)
  sub_data$risk_group <- ifelse(sub_data$risk_score > median_risk, "High Risk", "Low Risk")
  # 设置因子水平，保证 Low 在前，High 在后 (绘图时蓝色在前红色在后)
  sub_data$risk_group <- factor(sub_data$risk_group, levels = c("Low Risk", "High Risk"))
  
  # 如果分组失败(只有一组)，跳过
  if (length(unique(sub_data$risk_group)) < 2) return(NULL)
  
  # ==========================================================================
  # --- 步骤 3: 绘制 KM 曲线 (SCI 风格) ---
  # ==========================================================================
  km_filename <- paste0(cancer_name, "_", type_name, "_KM.pdf")
  km_full_path <- file.path(save_path, km_filename)
  
  # 构建公式
  km_form <- as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ risk_group"))
  fit_km <- survfit(km_form, data = sub_data)
  
  # !重要修复!: 手动注入 formula，防止 ggsurvplot 报错
  fit_km$call$formula <- km_form
  fit_km$call$data <- quote(sub_data)
  
  p_km <- ggsurvplot(
    fit_km,
    data = sub_data,
    
    # 统计参数
    pval = TRUE,             # 显示 P 值
    pval.method = TRUE,      # 显示检验方法 (Log-rank)
    conf.int = TRUE,         # 显示置信区间阴影
    conf.int.alpha = 0.1,    # 置信区间透明度
    
    # 风险表参数
    risk.table = TRUE,       # 显示风险表
    risk.table.col = "strata", # 风险表颜色跟随曲线
    risk.table.height = 0.25, # 风险表高度
    
    # 视觉风格
    palette = "jco",         # JCO 期刊经典配色 (蓝/黄)
    ggtheme = theme_classic() + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)),
    linetype = 1,            # 实线
    size = 1,                # 曲线粗细
    
    # 标签
    title = paste0(cancer_name, " - ", type_name, " Survival"),
    xlab = "Time (Days)",
    legend.title = "Risk Group",
    legend.labs = c("Low Risk", "High Risk")
  )
  
  # 保存 KM PDF
  pdf(km_full_path, width = 7, height = 7, onefile = FALSE)
  print(p_km, newpage = FALSE)
  dev.off()
  
  # ==========================================================================
  # --- 步骤 4: 绘制 ROC 曲线 (SCI 风格) ---
  # ==========================================================================
  roc_filename <- paste0(cancer_name, "_", type_name, "_ROC.pdf")
  roc_full_path <- file.path(save_path, roc_filename)
  
  max_time <- max(sub_data[[time_col]], na.rm = TRUE)
  times_eval <- c(365, 1095, 1825) # 1年, 3年, 5年
  # 对应 JCO/Lancet 常用三色 (金、蓝、红)
  colors <- c("#E7B800", "#2E9FDF", "#FC4E07") 
  
  # 开启 PDF 设备
  pdf(roc_full_path, width = 6, height = 6)
  
  # 设置边距
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # 初始化空画板
  plot(0, 0, type="n", xlim=c(0,1), ylim=c(0,1), 
       xlab="False Positive Rate (1 - Specificity)", 
       ylab="True Positive Rate (Sensitivity)",
       main=paste0(cancer_name, " - ", type_name, " ROC Curve"),
       cex.lab=1.2, cex.axis=1.2, font.lab=2)
  
  # 添加对角参考线
  abline(0, 1, col="gray80", lty=2, lwd=1.5)
  
  legend_labels <- c()
  
  # 循环计算 1, 3, 5 年
  for (i in 1:3) {
    t_point <- times_eval[i]
    if (t_point < max_time) {
      roc_res <- tryCatch({
        survivalROC(Stime = sub_data[[time_col]], 
                    status = sub_data[[status_col]], 
                    marker = sub_data$risk_score, 
                    predict.time = t_point, 
                    method = "NNE", 
                    span = 0.25 * nrow(sub_data)^(-0.20))
      }, error = function(e) NULL)
      
      if (!is.null(roc_res) && !is.na(roc_res$AUC)) {
        lines(roc_res$FP, roc_res$TP, col=colors[i], lwd=2.5) # 线条加粗
        auc_val <- sprintf("%.3f", roc_res$AUC) # 保留3位小数
        legend_labels <- c(legend_labels, paste0(i*1, "-Year AUC = ", auc_val))
      } else {
        legend_labels <- c(legend_labels, paste0(i*1, "-Year (NA)"))
      }
    } else {
      legend_labels <- c(legend_labels, paste0(i*1, "-Year (No Data)"))
    }
  }
  
  # 添加图例
  legend("bottomright", legend = legend_labels, col = colors, lwd = 2.5, 
         bty = "n", cex = 1.0, title = "AUC Value")
  
  dev.off() # 关闭设备
  
  return(TRUE)
}

# ==============================================================================
# 4. 主程序循环
# ==============================================================================

# 获取所有癌症类型并去除 NA
cancer_types <- sort(unique(TCGA_PANCAN_Merged_Final$CancerType))
cancer_types <- cancer_types[!is.na(cancer_types)]

cat("========================================================\n")
cat("开始批量绘图任务\n")
cat("结果将保存至:", file.path(getwd(), output_dir), "\n")
cat("========================================================\n")

count <- 0
total_cancers <- length(cancer_types)

for (cancer in cancer_types) {
  count <- count + 1
  cat(sprintf("[%02d/%02d] 正在处理 %s ... ", count, total_cancers, cancer))
  
  # 1. 为每个癌症创建单独的子文件夹
  cancer_dir <- file.path(output_dir, cancer)
  if (!dir.exists(cancer_dir)) dir.create(cancer_dir)
  
  # 2. 提取该癌症的数据
  cancer_data <- subset(TCGA_PANCAN_Merged_Final, CancerType == cancer)
  
  # 3. 循环 4 种生存数据类型
  plots_generated <- 0
  for (type_name in names(surv_types)) {
    cols <- surv_types[[type_name]]
    
    # 调用绘图函数
    res <- analyze_and_plot(
      data = cancer_data,
      cancer_name = cancer,
      type_name = type_name,
      time_col = cols[1],
      status_col = cols[2],
      gene_list = genes,
      save_path = cancer_dir
    )
    
    if (!is.null(res)) plots_generated <- plots_generated + 1
  }
  
  cat(paste0("完成 (生成 ", plots_generated * 2, " 张图)\n"))
}

cat("========================================================\n")
cat("全部任务结束！\n")