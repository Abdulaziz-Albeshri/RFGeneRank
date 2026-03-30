# ---- internal helpers (must appear BEFORE sign_importance) -------------------

.rfgr_as_matrix <- function(X, assay_name = NULL,
                            orientation = c("samples_by_genes","genes_by_samples")) {
  orientation <- match.arg(orientation)
  if (inherits(X, "SummarizedExperiment")) {
    aa <- SummarizedExperiment::assays(X)
    X  <- if (is.null(assay_name)) aa[[1]] else aa[[assay_name]]
  } else if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (orientation == "genes_by_samples") X <- t(X)
  storage.mode(X) <- "double"
  X
}

.rfgr_try_get_y <- function(fit) {
  y <- NULL
  try({ if (!is.null(fit$oof$y))  y <- fit$oof$y }, silent = TRUE)
  try({oof0 <- .rfgr_oof(fit)
  if (!is.null(oof0$y)) y <- oof0$y
  }, silent = TRUE)
  try({ if (!is.null(fit$trainingData$.outcome)) y <- fit$trainingData$.outcome }, silent = TRUE)
  try({ if (!is.null(fit$y)) y <- fit$y }, silent = TRUE)
  y
}

.rfgr_get_importance <- function(fit) {
  imp <- NULL

  # 0) RFGeneRank pipeline: importance is stored as a data.frame with columns gene, importance
  try({
    imp0 <- .rfgr_imp(fit)
    if (!is.null(imp0) && is.data.frame(imp0) &&
    all(c("gene","importance") %in% colnames(imp0))) {
    v <- imp0$importance
    names(v) <- imp0$gene
      imp <- v
    }
  }, silent = TRUE)

  # 1) ranger object stored in the model component (some wrappers)
  if (is.null(imp)) try({
    mdl1 <- .rfgr_model(fit)
    if (!is.null(mdl1$variable.importance))
    imp <- mdl1$variable.importance
  }, silent = TRUE)

  # 2) plain ranger object
  if (is.null(imp)) try({
    if (!is.null(fit$variable.importance))
      imp <- fit$variable.importance
  }, silent = TRUE)

  # 3) custom numeric vector
  if (is.null(imp)) try({
    if (!is.null(fit$importance) && is.numeric(fit$importance))
      imp <- fit$importance
  }, silent = TRUE)

  # return named numeric vector or NULL
  if (!is.null(imp)) {
  imp0 <- try(.rfgr_imp(fit), silent = TRUE)
  if (inherits(imp0, "try-error")) imp0 <- NULL

  if (is.null(names(imp)) && !is.null(imp0) && !is.null(imp0$gene)) {
    names(imp) <- imp0$gene
  }
  return(imp)
}
  NULL
}


.rfgr_safe_log2fc <- function(mean_case, mean_ctrl, eps = 1e-6) {
  mc <- as.numeric(mean_case); m0 <- as.numeric(mean_ctrl)
  log2((mc + eps) / (m0 + eps))
}

.rfgr_predict_proba_fun <- function(fit) {
  if (inherits(fit, "ranger")) {
    return(function(newdata) {
      p <- predict(fit, data = newdata, type = "response")$predictions
      if (is.matrix(p)) p[, 2, drop = TRUE] else as.numeric(p)
    })
  }
  if (inherits(fit, "randomForest")) {
    return(function(newdata) as.numeric(predict(fit, newdata, type = "prob")[, 2]))
  }
  if ("train" %in% class(fit)) {
    return(function(newdata) as.numeric(predict(fit, newdata, type = "prob")[, 2]))
  }
  if (!is.null(fit$predict_proba) && is.function(fit$predict_proba)) {
    return(function(newdata) as.numeric(fit$predict_proba(newdata)))
  }
  NULL
}

# ---- main API ---------------------------------------------------------------

#' Signed feature importance for RFGeneRank models
#'
#' @description 
#' Add direction (+/-) to RF importance using group means ("mean"),
#' external DE log2FC ("de"), or SHAP ("shap").
#'
#' @param fit Trained model object from rank_genes() or similar.
#' @param X   Expression (matrix/data.frame/SummarizedExperiment).
#' @param y   Factor labels (optional if retrievable from fit).
#' @param method c("mean","de","shap"). Default "mean".
#' @param de_table data.frame with columns gene, log2FC (for method="de").
#' @param case_level Positive class label (default = last level of y).
#' @param orientation "samples_by_genes" or "genes_by_samples".
#' @param assay_name SE assay name/index.
#' @param shap_n Integer; number of Monte Carlo samples for SHAP.
#' @param seed RNG seed.
#'
#' @return 
#' data.frame with gene, importance, direction, signed_importance,
#' mean_case, mean_ctrl, mean_diff, log2FC, shap_dir.
#' @importFrom stats relevel
#' @export
#' @examples
#' # Toy expression matrix: samples x genes
#' X <- matrix(
#'   rnorm(4 * 5),
#'   nrow = 4,
#'   dimnames = list(
#'     paste0("sample", 1:4),
#'     paste0("gene",   1:5)
#'   )
#' )
#'
#' # Binary labels
#' y <- factor(c("Control", "Control", "Case", "Case"))
#'
#' # Minimal "fit" object: a list with a numeric importance vector
#' fit <- list(
#'   importance = setNames(
#'     c(0.8, 0.6, 0.4, 0.2, 0.1),
#'     paste0("gene", 1:5)
#'   )
#' )
#'
#' res <- sign_importance(
#'   fit = fit,
#'   X   = X,
#'   y   = y,
#'   method = "mean"
#' )
#'
#' head(res)
#' # In practice, sign_importance() is called on a GeneRankFit object
#' # produced by the RFGeneRank workflow, for example:
#' #
#' #   fit <- gene_rank(se, genes = rownames(se), ...)
#' #   sig_imp <- sign_importance(fit, X = expr_matrix, y = outcome)
#' #   head(sig_imp)
sign_importance <- function(
  fit, X, y = NULL,
  method = c("mean","de","shap"),
  de_table = NULL,
  case_level = NULL,
  orientation = c("samples_by_genes","genes_by_samples"),
  assay_name = NULL,
  shap_n = 50,
  seed = 1
) {
  orientation <- match.arg(orientation)
  method      <- match.arg(method)

  # matrix + labels
  X_mat <- .rfgr_as_matrix(X, assay_name = assay_name, orientation = orientation)
  if (is.null(y)) y <- .rfgr_try_get_y(fit)
  if (is.null(y)) stop("`y` is required if it cannot be extracted from `fit`.")
  y <- as.factor(y)
  if (is.null(case_level)) case_level <- levels(y)[length(levels(y))]
  y  <- relevel(y, ref = setdiff(levels(y), case_level)[1])

  # align dims
  if (nrow(X_mat) != length(y)) {
    if (ncol(X_mat) == length(y)) {
      X_mat <- t(X_mat)
    } else stop("Dimensions of X and y do not align.")
  }

  # raw importance
  imp <- .rfgr_get_importance(fit)
  if (is.null(imp)) stop("Could not extract variable importance from `fit`.")
  imp <- imp[intersect(names(imp), colnames(X_mat))]
  imp <- imp[order(imp, decreasing = TRUE)]
  genes <- names(imp)
  Xg <- X_mat[, genes, drop = FALSE]

  dir_vec <- rep(NA_real_, length(genes)); names(dir_vec) <- genes

  if (method == "mean") {
    mean_case <- colMeans(Xg[y == case_level, , drop=FALSE], na.rm=TRUE)
    mean_ctrl <- colMeans(Xg[y != case_level, , drop=FALSE], na.rm=TRUE)
    mean_diff <- mean_case - mean_ctrl
    dir_vec   <- sign(mean_diff)
    log2FC    <- .rfgr_safe_log2fc(mean_case, mean_ctrl)
    shap_dir  <- rep(NA_real_, length(genes))
  } else if (method == "de") {
    if (is.null(de_table) || !all(c("gene","log2FC") %in% colnames(de_table)))
      stop("For method='de', provide de_table with columns: gene, log2FC")
    de_map  <- setNames(de_table$log2FC, de_table$gene)
    log2FC  <- de_map[genes]; log2FC[is.na(log2FC)] <- 0
    dir_vec <- sign(log2FC)
    mean_case <- colMeans(Xg[y == case_level, , drop=FALSE], na.rm=TRUE)
    mean_ctrl <- colMeans(Xg[y != case_level, , drop=FALSE], na.rm=TRUE)
    mean_diff <- mean_case - mean_ctrl
    shap_dir  <- rep(NA_real_, length(genes))
  } else { # shap
    if (!requireNamespace("fastshap", quietly = TRUE))
      stop("method='shap' requires the 'fastshap' package. Install with install.packages('fastshap').")
    pred_fun <- .rfgr_predict_proba_fun(fit)
    if (is.null(pred_fun)) stop("Could not build probability predictor for SHAP.")

    S <- fastshap::explain(
      object = fit,
      X = as.data.frame(Xg),
      pred_wrapper = function(object, newdata) pred_fun(newdata),
      nsim = shap_n,
      adjust = TRUE
    )
    shap_dir  <- sign(colMeans(as.matrix(S), na.rm = TRUE))
    dir_vec   <- shap_dir
    mean_case <- colMeans(Xg[y == case_level, , drop=FALSE], na.rm=TRUE)
    mean_ctrl <- colMeans(Xg[y != case_level, , drop=FALSE], na.rm=TRUE)
    mean_diff <- mean_case - mean_ctrl
    log2FC    <- .rfgr_safe_log2fc(mean_case, mean_ctrl)
  }

  dir_vec[is.na(dir_vec)] <- 0
  signed_imp <- as.numeric(imp) * dir_vec

  out <- data.frame(
    gene = genes,
    importance = as.numeric(imp),
    direction  = as.integer(dir_vec),
    signed_importance = signed_imp,
    mean_case = mean_case[genes],
    mean_ctrl = mean_ctrl[genes],
    mean_diff = mean_diff[genes],
    log2FC    = log2FC,
    shap_dir  = shap_dir[genes],
    row.names = NULL,
    check.names = FALSE
  )
  out[order(-abs(out$signed_importance)), ]
}
