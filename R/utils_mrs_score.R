# Utility functions for the final 10-gene MRS score.

source(file.path(getwd(), "R", "config.R"))

calculate_mrs_score <- function(expr, coefficients = mrs_coefficients,
                                gene_as = c("row", "column"),
                                missing = c("error", "zero")) {
  gene_as <- match.arg(gene_as)
  missing <- match.arg(missing)

  if (gene_as == "column") {
    expr_use <- as.matrix(expr[, names(coefficients), drop = FALSE])
  } else {
    expr_use <- t(as.matrix(expr[names(coefficients), , drop = FALSE]))
  }

  missing_genes <- setdiff(names(coefficients), colnames(expr_use))
  if (length(missing_genes) > 0) {
    if (missing == "error") {
      stop("Missing genes in expression matrix: ", paste(missing_genes, collapse = ", "))
    }
    for (g in missing_genes) expr_use[, g] <- 0
  }

  expr_use <- expr_use[, names(coefficients), drop = FALSE]
  as.numeric(expr_use %*% coefficients)
}

median_risk_group <- function(score) {
  ifelse(score > stats::median(score, na.rm = TRUE), "High", "Low")
}
