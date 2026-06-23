# Curated release script for: MRS-LUAD
# Source script: cox/tcgaluad10gene.R
# Role: Summarize the 10-gene model in TCGA-LUAD.
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
# 单个队列：10个基因的单基因KM曲线，2行5列
# 只需修改 DATA_OBJ_NAME
# ==============================================================================

library(survival)
library(survminer)
library(ggplot2)
library(cowplot)

# ==============================================================================
# 1. 参数设置
# ==============================================================================
DATA_OBJ_NAME <- "tcgaluad"   # 改成你的数据对象名，例如 "GSE31210"、"GSE72094"、"tcgaluad"、"tcgalusc"

target_genes <- c(
  "RRAS", "GPM6B", "FMO3", "CHRDL1", "ADRB2",
  "LPL", "ANGPT1", "SPOCK2", "P2RY14", "ADH1B"
)

km_colors <- c("High expression" = "#E64B35", "Low expression" = "#4DBBD5")

out_pdf <- paste0("KM_2x5_", DATA_OBJ_NAME, "_10genes.pdf")
out_csv <- paste0("KM_2x5_", DATA_OBJ_NAME, "_10genes_stats.csv")

# ==============================================================================
# 2. 安全读取数据对象
# ==============================================================================
if (!exists(DATA_OBJ_NAME)) {
  stop(paste0("对象不存在：", DATA_OBJ_NAME))
}
dat0 <- get(DATA_OBJ_NAME)

cat("Dataset object name: ", DATA_OBJ_NAME, "\n")
cat("Dimension:\n")
print(dim(dat0))
cat("Columns preview:\n")
print(head(colnames(dat0), 30))

# ==============================================================================
# 3. 检查生存列
# ==============================================================================
required_surv_cols <- c("OS.time", "OS")
if (!all(required_surv_cols %in% colnames(dat0))) {
  stop("数据中缺少 OS.time 或 OS 列。")
}

# ==============================================================================
# 4. 单基因KM函数
# ==============================================================================
make_single_gene_km <- function(data, gene, palette) {
  
  empty_plot <- ggplot() +
    theme_void() +
    ggtitle(paste0(gene, "\nNo valid data"))
  
  # 检查基因是否存在
  if (!(gene %in% colnames(data))) {
    return(list(
      plot = empty_plot,
      stats = data.frame(
        Gene = gene,
        N = NA,
        Cutoff = NA,
        High_n = NA,
        Low_n = NA,
        HR = NA,
        HR_lower95 = NA,
        HR_upper95 = NA,
        Cox_P = NA,
        Logrank_P = NA,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # 提取所需列并去NA
  df <- data[, c("OS.time", "OS", gene)]
  df <- na.omit(df)
  
  # 样本太少就跳过
  if (nrow(df) < 10) {
    return(list(
      plot = empty_plot,
      stats = data.frame(
        Gene = gene,
        N = nrow(df),
        Cutoff = NA,
        High_n = NA,
        Low_n = NA,
        HR = NA,
        HR_lower95 = NA,
        HR_upper95 = NA,
        Cox_P = NA,
        Logrank_P = NA,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # 表达无变异也跳过
  if (length(unique(df[[gene]])) < 2) {
    return(list(
      plot = empty_plot,
      stats = data.frame(
        Gene = gene,
        N = nrow(df),
        Cutoff = NA,
        High_n = NA,
        Low_n = NA,
        HR = NA,
        HR_lower95 = NA,
        HR_upper95 = NA,
        Cox_P = NA,
        Logrank_P = NA,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # 按中位数分组
  cutoff <- median(df[[gene]], na.rm = TRUE)
  df$expr_group <- ifelse(df[[gene]] > cutoff, "High expression", "Low expression")
  df$expr_group <- factor(df$expr_group, levels = c("High expression", "Low expression"))
  
  high_n <- sum(df$expr_group == "High expression")
  low_n  <- sum(df$expr_group == "Low expression")
  
  # 如果某组为空，则跳过
  if (high_n == 0 || low_n == 0) {
    return(list(
      plot = empty_plot,
      stats = data.frame(
        Gene = gene,
        N = nrow(df),
        Cutoff = cutoff,
        High_n = high_n,
        Low_n = low_n,
        HR = NA,
        HR_lower95 = NA,
        HR_upper95 = NA,
        Cox_P = NA,
        Logrank_P = NA,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # KM 拟合
  fit <- tryCatch(
    survfit(Surv(OS.time, OS) ~ expr_group, data = df),
    error = function(e) NULL
  )
  
  # Cox（给HR）
  cox_fit <- tryCatch(
    coxph(Surv(OS.time, OS) ~ expr_group, data = df),
    error = function(e) NULL
  )
  
  hr <- NA
  hr_l <- NA
  hr_u <- NA
  cox_p <- NA
  
  if (!is.null(cox_fit)) {
    sm <- summary(cox_fit)
    if (nrow(sm$coefficients) >= 1) {
      hr    <- sm$coefficients[1, "exp(coef)"]
      cox_p <- sm$coefficients[1, "Pr(>|z|)"]
      hr_l  <- sm$conf.int[1, "lower .95"]
      hr_u  <- sm$conf.int[1, "upper .95"]
    }
  }
  
  # log-rank
  logrank_p <- NA
  lr <- tryCatch(
    survdiff(Surv(OS.time, OS) ~ expr_group, data = df),
    error = function(e) NULL
  )
  if (!is.null(lr)) {
    logrank_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  }
  
  # 画图
  if (is.null(fit)) {
    p <- empty_plot
  } else {
    p_obj <- ggsurvplot(
      fit,
      data = df,
      pval = TRUE,
      conf.int = FALSE,
      risk.table = FALSE,
      palette = palette,
      legend.title = gene,
      legend.labs = c("High expression", "Low expression"),
      title = gene,
      xlab = "Time",
      ggtheme = theme_bw() +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          legend.position = "top",
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8)
        )
    )
    
    p <- p_obj$plot +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = paste0(
          "HR=", ifelse(is.na(hr), "NA", sprintf("%.2f", hr)),
          "\n95%CI ", ifelse(is.na(hr_l), "NA", sprintf("%.2f", hr_l)),
          "-", ifelse(is.na(hr_u), "NA", sprintf("%.2f", hr_u))
        ),
        hjust = 1.05, vjust = 1.5,
        size = 3.0,
        fontface = "bold"
      )
  }
  
  stats_df <- data.frame(
    Gene = gene,
    N = nrow(df),
    Cutoff = cutoff,
    High_n = high_n,
    Low_n = low_n,
    HR = hr,
    HR_lower95 = hr_l,
    HR_upper95 = hr_u,
    Cox_P = cox_p,
    Logrank_P = logrank_p,
    stringsAsFactors = FALSE
  )
  
  return(list(plot = p, stats = stats_df))
}

# ==============================================================================
# 5. 批量生成10个基因图
# ==============================================================================
plot_list <- list()
stats_list <- list()

cat("Generating KM plots...\n")

for (g in target_genes) {
  cat("Processing: ", g, "\n")
  res <- make_single_gene_km(dat0, g, km_colors)
  plot_list[[g]] <- res$plot
  stats_list[[g]] <- res$stats
}

stats_df <- do.call(rbind, stats_list)

cat("\nKM stats summary:\n")
print(stats_df)

write.csv(stats_df, out_csv, row.names = FALSE)

# ==============================================================================
# 6. 拼图：2行5列
# ==============================================================================
final_plot <- plot_grid(
  plotlist = plot_list,
  ncol = 5,
  nrow = 2,
  align = "hv"
)

ggsave(
  filename = out_pdf,
  plot = final_plot,
  width = 22,
  height = 9,
  limitsize = FALSE
)

cat("\nDone!\n")
cat("PDF saved to: ", out_pdf, "\n")
cat("Stats saved to: ", out_csv, "\n")