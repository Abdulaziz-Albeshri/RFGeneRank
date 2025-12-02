# ---- R/globals.R ------------------------------------------------------------

# CRAN NOTE fix when using non-standard evaluation in ggplot2/dplyr, etc.
utils::globalVariables(c(
  ".", "gene", "importance", "SelectedInFolds", "freq", "score", "label",
  "UMAP1", "UMAP2", "PC1", "PC2", "state", "x", "y", "xend", "yend",
  "signed_importance", "direction", "fdr", "pct", "Freq", "True", "Predicted", "method", "group"
))

# Also declare fpr/tpr for plot_roc()
utils::globalVariables(c("fpr","tpr","FPR","TPR","SHAP"))

#' @importFrom stats addmargins na.omit reformulate setNames model.matrix as.formula median quantile relevel
#' @importFrom methods new
#' @importFrom stats addmargins na.omit reformulate setNames model.matrix as.formula
#' @importFrom stats median quantile relevel lm anova cor.test resid loess p.adjust
#' @importFrom utils tail
NULL
