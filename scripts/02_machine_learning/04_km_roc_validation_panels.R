# Curated release script for: MRS-LUAD
# Source script: intersection_mechine/mechine/outside_data/KM曲线发表文章重做.R
# Role: Generate Kaplan-Meier and time-dependent ROC panels for validation cohorts.
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
# 1. 环境设置与库加载
# ==============================================================================
library(survival)
library(survminer)
library(ggplot2)
library(survivalROC)
library(cowplot)  # 用于拼图
library(grid)     # 用于绘制空白占位符

# ==============================================================================
# 2. 参数定义
# ==============================================================================
# 基因列表
target_genes <- c("RRAS", "GPM6B", "FMO3", "CHRDL1", "ADRB2", 
                  "LPL", "ANGPT1", "SPOCK2", "P2RY14", "ADH1B")

# KM 配色方案
risk_colors <- c("#E64B35", "#4DBBD5")

# ROC 配色方案 (用于曲线和文字)
# 确保名称与后面代码中的 label 一致
roc_curve_colors <- c("1-Year" = "#E64B35", "3-Year" = "#00A087", "5-Year" = "#3C5488")

# ==============================================================================
# 3. 核心函数：生成 KM 和 ROC 对象
# ==============================================================================
get_plot_objects <- function(data, dataset_name) {
  
  # --- 空数据处理：返回空白图 ---
  empty_plot <- ggplot() + theme_void()
  res_list <- list(km = empty_plot, roc = empty_plot)
  
  if (is.null(data)) return(res_list)
  
  # --- 数据清洗 ---
  required_cols <- c("OS.time", "OS", target_genes)
  if (!all(required_cols %in% colnames(data))) return(res_list)
  data_clean <- na.omit(data[, required_cols])
  if (nrow(data_clean) < 10) return(res_list)
  
  # --- Cox 模型与风险评分 ---
  formula_str <- paste("Surv(OS.time, OS) ~", paste(target_genes, collapse = " + "))
  cox_model <- tryCatch({
    coxph(as.formula(formula_str), data = data_clean, control = coxph.control(iter.max = 100))
  }, error = function(e) return(NULL))
  
  if (is.null(cox_model) || any(is.na(coef(cox_model)))) return(res_list)
  
  # 计算 Risk Score
  coefficients <- coef(cox_model)
  expr_matrix <- as.matrix(data_clean[, target_genes])
  data_clean$risk_score <- as.vector(expr_matrix %*% coefficients)
  
  # 分组
  median_risk <- median(data_clean$risk_score, na.rm = TRUE)
  data_clean$risk_group <- ifelse(data_clean$risk_score > median_risk, "High Risk", "Low Risk")
  data_clean$risk_group <- factor(data_clean$risk_group, levels = c("High Risk", "Low Risk"))
  
  # --- 生成 KM Plot 对象 ---
  fit_risk <- survfit(Surv(OS.time, OS) ~ risk_group, data = data_clean)
  
  km_obj <- ggsurvplot(
    fit_risk, 
    data = data_clean, 
    pval = TRUE, 
    conf.int = TRUE, 
    risk.table = FALSE, 
    palette = risk_colors,
    title = dataset_name,
    xlab = "Time(days)",
    ggtheme = theme_bw() + theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = "top"
    )
  )
  
  # --- 生成 ROC Plot 对象 (修改了这里) ---
  max_time <- max(data_clean$OS.time, na.rm = TRUE)
  target_times <- c(365, 1095, 1825)
  check_times <- pmin(target_times, max_time - 1)
  time_labels <- c("1-Year", "3-Year", "5-Year")
  
  roc_data_list <- list()
  auc_map <- list() # 用于存储AUC值，方便后续按颜色提取
  
  for (i in 1:3) {
    t <- check_times[i]
    lbl <- time_labels[i]
    roc_res <- tryCatch({
      survivalROC(Stime = data_clean$OS.time, status = data_clean$OS, 
                  marker = data_clean$risk_score, predict.time = t, 
                  method = "NNE", span = 0.25 * nrow(data_clean)^(-0.20))
    }, error = function(e) return(NULL))
    
    if (!is.null(roc_res) && !is.na(roc_res$AUC)) {
      roc_data_list[[lbl]] <- data.frame(FP = roc_res$FP, TP = roc_res$TP, Time = lbl)
      # 存储格式化的 AUC 字符串
      auc_map[[lbl]] <- sprintf("%.2f", roc_res$AUC)
    }
  }
  
  if (length(roc_data_list) > 0) {
    roc_df_all <- do.call(rbind, roc_data_list)
    roc_df_all$Time <- factor(roc_df_all$Time, levels = time_labels)
    
    # 基础绘图
    roc_gg <- ggplot(roc_df_all, aes(x = FP, y = TP, color = Time)) +
      geom_path(size = 1) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
      scale_color_manual(values = roc_curve_colors) +
      labs(x = "1 - Specificity", y = "Sensitivity", title = paste0("ROC: ", dataset_name)) +
      theme_bw() +
      theme(
        panel.grid = element_blank(), 
        legend.position = "none", # 隐藏默认图例
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
      )
    
    # 【核心修改】：手动循环添加彩色文字注释
    # 初始 Y 轴位置 (从 0.25 开始向下排列，根据你的图片位置调整)
    current_y <- 0.25 
    y_gap <- 0.08    # 行间距
    
    # 按 1 -> 3 -> 5 的顺序添加
    for (lbl in time_labels) {
      if (!is.null(auc_map[[lbl]])) {
        # 构造文本: e.g. "1-Year: 0.83"
        label_text <- paste0(lbl, ": ", auc_map[[lbl]])
        
        # 获取对应颜色
        txt_color <- roc_curve_colors[[lbl]]
        
        # 添加文本层
        roc_gg <- roc_gg + 
          annotate("text", x = 0.98, y = current_y, 
                   label = label_text, 
                   color = txt_color, 
                   hjust = 1, # 右对齐
                   size = 3.5, 
                   fontface = "bold")
        
        # 更新 Y 坐标，让下一行往下移
        current_y <- current_y - y_gap
      }
    }
    
  } else {
    roc_gg <- empty_plot
  }
  
  return(list(km = km_obj$plot, roc = roc_gg))
}

# ==============================================================================
# 4. 准备数据列表（更新了命名）
# ==============================================================================

# 定义一个辅助函数来安全获取数据，如果不存在返回NULL
safe_get <- function(var_name) {
  if (exists(var_name)) return(get(var_name)) else return(NULL)
}

# --- 修改了名字：Yes Smoke -> Smoke ---

# 第一组：GSE31210
row1_names <- c("GSE31210 (Total)", "GSE31210 (No Smoke)", "GSE31210 (Smoke)")
row1_data  <- list(
  safe_get("GSE31210"), 
  safe_get("GSE31210_smoke_no"), 
  safe_get("GSE31210_smoke_yes")
)

# 第二组：GSE72094
row3_names <- c("GSE72094 (Total)", "GSE72094 (No Smoke)", "GSE72094 (Smoke)")
row3_data  <- list(
  safe_get("GSE72094"), 
  safe_get("GSE72094_smoke_no"), 
  safe_get("GSE72094_smoke_yes")
)

# 第三组：TCGA
row5_names <- c("TCGA-LUAD", "TCGA-LUSC", "")
row5_data  <- list(
  safe_get("tcgaluad"), 
  safe_get("tcgalusc"), 
  NULL # 占位符
)

# ==============================================================================
# 5. 生成所有图对象 (保持不变)
# ==============================================================================

cat("正在生成图表对象...\n")

# 容器
plots_r1_km <- list(); plots_r2_roc <- list()
plots_r3_km <- list(); plots_r4_roc <- list()
plots_r5_km <- list(); plots_r6_roc <- list()

# 处理第一行/第二行 (GSE31210)
for(i in 1:3) {
  res <- get_plot_objects(row1_data[[i]], row1_names[i])
  plots_r1_km[[i]] <- res$km
  plots_r2_roc[[i]] <- res$roc
}

# 处理第三行/第四行 (GSE72094)
for(i in 1:3) {
  res <- get_plot_objects(row3_data[[i]], row3_names[i])
  plots_r3_km[[i]] <- res$km
  plots_r4_roc[[i]] <- res$roc
}

# 处理第五行/第六行 (TCGA)
for(i in 1:3) {
  res <- get_plot_objects(row5_data[[i]], row5_names[i])
  plots_r5_km[[i]] <- res$km
  plots_r6_roc[[i]] <- res$roc
}

# ==============================================================================
# 6. 拼图并输出高清PDF
# ==============================================================================

cat("正在排版并输出 PDF...\n")

# Row 1: GSE31210 KM
g1 <- plot_grid(plotlist = plots_r1_km, ncol = 3, align = "h")
# Row 2: GSE31210 ROC
g2 <- plot_grid(plotlist = plots_r2_roc, ncol = 3, align = "h")
# Row 3: GSE72094 KM
g3 <- plot_grid(plotlist = plots_r3_km, ncol = 3, align = "h")
# Row 4: GSE72094 ROC
g4 <- plot_grid(plotlist = plots_r4_roc, ncol = 3, align = "h")
# Row 5: TCGA KM
g5 <- plot_grid(plotlist = plots_r5_km, ncol = 3, align = "h")
# Row 6: TCGA ROC
g6 <- plot_grid(plotlist = plots_r6_roc, ncol = 3, align = "h")

# 将6行组合在一起
# 相对高度设置：KM图 (1.2) 比 ROC图 (1.0) 稍高一点，视觉上更平衡
final_grid <- plot_grid(g1, g2, g3, g4, g5, g6, 
                        ncol = 1, 
                        rel_heights = c(1.2, 1, 1.2, 1, 1.2, 1))

# 输出
pdf_filename <- "Survival_Matrix_Report_Fixed.pdf"
ggsave(pdf_filename, final_grid, width = 15, height = 22, limitsize = FALSE)

cat(paste0("完成！\n1. AUC颜色已修复 \n2. 标题已更新(Yes Smoke -> Smoke)\n文件已保存至: ", pdf_filename, "\n"))
