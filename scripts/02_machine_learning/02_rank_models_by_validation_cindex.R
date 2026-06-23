# Curated release script for: MRS-LUAD
# Source script: intersection_mechine/mechine/outside_data/plot_cindex_validate_only.R
# Role: Rank models by validation-set C-index and draw the model-performance heatmap.
#
# Usage note:
#   Run scripts from the repository root when possible. Public datasets and large
#   intermediate R objects are not bundled. Set MRS_PROJECT_ROOT and SCFEA_ROOT
#   to local analysis directories before rerunning on another system.

if (!exists("mrs_project_root")) {
  config_file <- file.path(getwd(), "R", "config.R")
  if (file.exists(config_file)) source(config_file)
}

root <- "~/lung_cancer/intersection_mechine/mechine/outside_data"
setwd(path.expand(root))

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(aplot))

load("workspace.RData")

if (!exists("res55")) {
  stop("Object res55 was not found in workspace.RData")
}
if (!"Cindex.res" %in% names(res55)) {
  stop("res55 does not contain Cindex.res")
}

validate_set <- c("TCGA_LUAD_Vali", "GSE31210_Vali", "GSE72094_Vali")
plot_order <- validate_set

cindex_dis_validate_only <- function(object, color = NULL, dataset_col = NULL,
                                     validate_set, order = NULL,
                                     width = NULL, height = NULL) {
  if (is.null(width) == TRUE) {
    width <- 0.35
  } else {
    width <- width
  }
  if (is.null(height) == TRUE) {
    height <- 0.01
  } else {
    height <- height
  }
  if (is.null(color) == TRUE) {
    color <- c("#0084A7", "#F5FACD", "#E05D00", "#79AF97", "#8491B4")
  } else {
    color <- color
  }
  if (is.null(dataset_col) == TRUE) {
    dataset_col <- c("#3182BDFF", "#E6550DFF", "#31A354FF",
                     "#756BB1FF", "#636363FF", "#6BAED6FF",
                     "#FD8D3CFF", "#74C476FF", "#9E9AC8FF",
                     "#969696FF", "#9ECAE1FF", "#FDAE6BFF",
                     "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF",
                     "#C6DBEFFF", "#FDD0A2FF", "#C7E9C0FF",
                     "#DADAEBFF", "#D9D9D9FF")
  } else {
    dataset_col <- dataset_col
  }

  cindex <- object[["Cindex.res"]]
  cindex$Cindex <- as.numeric(sprintf("%.2f", cindex$Cindex))

  cindex_vaidate <- cindex[cindex$ID %in% validate_set, ]
  mean_cindex_validate <- aggregate(x = cindex_vaidate$Cindex,
                                    by = list(cindex_vaidate$Model), FUN = mean)
  colnames(mean_cindex_validate) <- c("Model", "mean")
  mean_cindex_validate$mean <- as.numeric(sprintf("%.3f", mean_cindex_validate$mean))
  mean_cindex_validate <- mean_cindex_validate[order(mean_cindex_validate$mean,
                                                     decreasing = FALSE), ]
  mean_cindex_validate$Value <- "Mean C-index in validate cohorts"

  labels <- as.data.frame(unique(cindex$ID))
  colnames(labels) <- "Cohort"
  if (is.null(order) == TRUE) {
    labels <- labels
  } else {
    order <- order
    labels$Cohort <- factor(labels$Cohort, levels = order)
    cindex$ID <- factor(cindex$ID, levels = order)
  }

  cindex$Model <- factor(cindex$Model, levels = mean_cindex_validate$Model)
  mean_cindex_validate$Model <- factor(mean_cindex_validate$Model,
                                       levels = mean_cindex_validate$Model)

  p1 <- ggplot(cindex, aes(x = ID, y = Model)) +
    geom_tile(aes(fill = Cindex), color = "white", size = 0.5, alpha = 0.80) +
    geom_text(aes(label = Cindex), hjust = 0.5, vjust = 0.5, color = "black", size = 2.5) +
    scale_fill_gradient2(low = color[1], mid = color[2], high = color[3],
                         midpoint = median(cindex$Cindex), name = "C-index",
                         limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_text(size = 9),
          panel.grid = element_blank(), panel.background = element_rect(fill = "white"),
          axis.text.x = element_blank())

  p2 <- ggplot(labels, aes(Cohort, y = 1)) +
    geom_tile(aes(fill = Cohort), color = "white", size = 0.5) +
    scale_fill_manual(values = dataset_col, name = "Cohort") +
    theme_void()

  p4 <- ggplot(data = mean_cindex_validate, aes(x = Model, y = mean, fill = Value)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = color[5], name = "") +
    geom_text(aes(y = mean / 2, label = mean),
              vjust = 0.5, hjust = 0.5, size = 2.5) +
    theme(axis.title = element_text(size = 8),
          axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank(), panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "white")) +
    labs(y = "") +
    coord_flip()

  print(p1 %>% insert_top(p2, height = height) %>% insert_right(p4, width = width))
}

available_sets <- unique(as.character(res55[["Cindex.res"]]$ID))
cat("Available C-index cohorts:\n")
print(available_sets)
cat("Validate cohorts used for mean C-index:\n")
print(validate_set)

pdf_file <- "cindex_validate_only_res55.pdf"
png_file <- "cindex_validate_only_res55.png"

if (capabilities("cairo")) {
  cairo_pdf(pdf_file, width = 9, height = 18)
} else {
  pdf(pdf_file, width = 9, height = 18, useDingbats = FALSE)
}
cindex_dis_validate_only(res55, validate_set = validate_set, order = plot_order,
                         width = 0.35, height = 0.01)
dev.off()

png(png_file, width = 2700, height = 5400, res = 300, type = "cairo")
cindex_dis_validate_only(res55, validate_set = validate_set, order = plot_order,
                         width = 0.35, height = 0.01)
dev.off()

cat("Saved files:\n")
print(normalizePath(c(pdf_file, png_file), mustWork = FALSE))