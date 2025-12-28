# ---- R/rank_genes.R ---------------------------------------------
# internal helpers (not exported)
#' @keywords internal
#' @noRd
.build_folds <- function(meta, label_col, batch_col = NULL, k = 5,
                         cv = c("kfold","lobo","groupk")) {
  cv <- match.arg(cv)
  n <- nrow(meta); idx <- seq_len(n)

  if (cv == "kfold") {
    y <- as.factor(meta[[label_col]])
    if (requireNamespace("caret", quietly = TRUE)) {
      tr_lists <- caret::createFolds(y, k = k, returnTrain = TRUE)
      return(lapply(seq_along(tr_lists), function(i){
        tr <- tr_lists[[i]]; te <- setdiff(idx, tr)
        list(train_idx = tr, test_idx = te, tag = paste0("fold=", i))
      }))
    } else {
      s <- sample(idx)
      bins <- split(s, cut(seq_along(s), breaks = k, labels = FALSE))
      return(lapply(seq_along(bins), function(i){
        te <- bins[[i]]; tr <- setdiff(idx, te)
        list(train_idx = tr, test_idx = te, tag = paste0("fold=", i))
      }))
    }
  }

  if (is.null(batch_col) || !(batch_col %in% names(meta)))
    stop("For cv='", cv, "' you must supply a valid batch_col in metadata.")

  grp <- as.factor(meta[[batch_col]]); ub <- levels(grp)

  if (cv == "lobo") {
    return(lapply(ub, function(b){
      te <- which(grp == b); tr <- setdiff(idx, te)
      list(train_idx = tr, test_idx = te, tag = paste0("test=", b))
    }))
  }

  # groupk
  k <- min(k, length(ub))
  if (requireNamespace("caret", quietly = TRUE)) {
    tr_lists <- caret::groupKFold(group = grp, k = k)
    return(lapply(seq_along(tr_lists), function(i){
      tr <- tr_lists[[i]]; te <- setdiff(idx, tr)
      list(train_idx = tr, test_idx = te, tag = paste0("fold=", i))
    }))
  } else {
    s <- sample(ub)
    chunks <- split(s, cut(seq_along(s), breaks = k, labels = FALSE))
    lapply(seq_along(chunks), function(i){
      te_g <- chunks[[i]]
      te <- which(grp %in% te_g); tr <- setdiff(idx, te)
      list(train_idx = tr, test_idx = te, tag = paste0("fold=", i))
    })
  }
}

#' @keywords internal
#' @noRd
.drop_constant_train <- function(X_tr, X_te, eps = 1e-12) {
  keep <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowVars(X_tr) > eps
  } else {
    apply(X_tr, 1L, var) > eps
  }
  list(X_tr = X_tr[keep, , drop = FALSE],
       X_te = X_te[keep, , drop = FALSE])
}

# -----------------------------------------------------------------------------#
#' Rank genes with batch-aware cross-validation (UMAP-free)
#'
#' Performs CV on a SummarizedExperiment, aggregates fold-normalized feature
#' importances, and stores out-of-fold predictions. Supports K-fold, LOBO
#' (Leave-One-Batch-Out), and group-K by batch. Batch correction, filtering,
#' transforms, and standardization are applied *inside folds* using TRAIN-only
#' statistics. If `auto_confounds=TRUE`, the function inspects batch~label
#' association and will switch to LOBO and enable fold-safe batch correction
#' if confounding is moderate/severe.
#'
#' @param se SummarizedExperiment
#' @param label_col character; class-label column in colData(se)
#' @param n_top integer; per-fold top-variance feature count (0 = off)
#' @param k integer; number of folds (ignored by LOBO)
#' @param trees integer; number of trees for ranger
#' @param importance "permutation" or "impurity"
#' @param class_weights named numeric vector or NULL (auto-computed if NULL and imbalance >= 1.5x)
#' @param threads integer; CPU threads for ranger
#' @param seed integer; RNG seed
#' @param fold_batch_correction logical; if TRUE, do train-only batch removal per fold
#' @param batch_col optional; name of batch column in colData(se)
#' @param batch_covariates optional character vector of covariates for batch model
#' @param filter_low_expr logical; drop genes expressed (>0) in < min_prop of TRAIN
#' @param min_prop numeric in (0,1]; minimum TRAIN proportion to keep a gene
#' @param transform "none" or "log1p"
#' @param standardize logical; z-score by TRAIN mean/SD (applied to train+test)
#' @param cv "kfold", "lobo", or "groupk"
#' @param auto_confounds logical; if TRUE, auto-switch to LOBO and enable fold batch-correction when confounded
#' @return GeneRankFit S4 object
#' @export
#' @examples
#' # For reproducibility, specify a fixed seed (e.g., set.seed(1)) before running this example.
#'
#' # Toy expression: 20 genes × 12 samples
#' expr <- matrix(stats::rnorm(20 * 12), nrow = 20)
#' rownames(expr) <- paste0("gene", 1:20)
#' colnames(expr) <- paste0("sample", 1:12)
#'
#' # Binary labels
#' label <- factor(rep(c("A", "B"), each = 6))
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(expr = expr),
#'   colData = data.frame(label = label)
#' )
#'
#' # Rank genes with default settings
#' rg <- rank_genes(
#'   se        = se,
#'   label_col = "label"
#' )
#'
#' str(rg)
rank_genes <- function(
  se, label_col = "state",
  n_top = 500, k = 5, trees = 1000,
  importance = c("permutation","impurity"),
  class_weights = NULL,
  threads = max(1, parallel::detectCores() - 1),
  seed = 1,
  fold_batch_correction = FALSE,
  batch_col = NULL,
  batch_covariates = NULL,
  filter_low_expr = FALSE,
  min_prop = 0.20,
  transform = c("none","log1p"),
  standardize = FALSE,
  cv = c("kfold","lobo","groupk"),
  auto_confounds = TRUE
) {
  importance <- match.arg(importance)
  transform  <- match.arg(transform)
  cv         <- match.arg(cv)

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("SummarizedExperiment is required.", call. = FALSE)
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("ranger is required.", call. = FALSE)

  expr <- SummarizedExperiment::assay(se, "expr")
meta <- as.data.frame(SummarizedExperiment::colData(se))
if (!label_col %in% names(meta)) stop("label_col not found in colData(se).", call. = FALSE)
y <- factor(meta[[label_col]])
if (nlevels(y) < 2L) stop("Need at least 2 classes.", call. = FALSE)
names(y) <- colnames(expr)

# Prevent double batch correction: if prepare_data() already corrected globally,
# do not allow in-fold correction here.
pd <- S4Vectors::metadata(se); if (is.null(pd)) pd <- list()
pd_batch_done <- isTRUE(pd$`..batch_corrected`) ||
                 identical(pd$`..batch_correction_scope`, "global")

if (pd_batch_done && isTRUE(fold_batch_correction)) {
  stop("Batch correction appears to have been applied in prepare_data(); ",
       "do NOT enable in-fold correction in rank_genes(). Pick exactly one path.")
}


# If prepare_data saved defaults, recover them (optional) ...


# If prepare_data saved defaults, recover them (optional)
if (is.null(batch_col) && "..batch_col" %in% names(meta)) {
  bc <- unique(na.omit(meta$..batch_col))
  batch_col <- if (length(bc) == 1 && !is.na(bc)) bc else NULL
}
if (is.null(batch_covariates) && "..batch_covariates" %in% names(meta)) {
  first <- meta$..batch_covariates[[1]]
  if (!is.null(first)) batch_covariates <- first
}

  # Confounding guardrail (internal helper)
  if (!is.null(batch_col) && batch_col %in% names(meta)) {
    conf <- .check_confounding(meta[[label_col]], meta[[batch_col]])
    if (!is.null(conf$msg)) message(conf$msg)
    if (isTRUE(auto_confounds) && isTRUE(conf$warn)) {
      if (cv != "lobo") {
        message("Auto-confounds: switching cv='lobo'.")
        cv <- "lobo"
      }
      if (!isTRUE(fold_batch_correction)) {
        message("Auto-confounds: enabling fold_batch_correction=TRUE.")
        fold_batch_correction <- TRUE
      }
    }
  }

  # Auto class weights (internal helper)
  if (is.null(class_weights)) {
    tab <- table(y)
    if (length(tab) >= 2 && (max(tab) / min(tab) >= 1.5)) {
      class_weights <- .compute_class_weights(y)
    }
  }

  # Build folds
  folds <- .build_folds(meta, label_col = label_col, batch_col = batch_col, k = k, cv = cv)

  sample_ids <- colnames(expr)
  oof_prob <- matrix(NA_real_, nrow = length(sample_ids), ncol = nlevels(y),
                     dimnames = list(sample_ids, levels(y)))
  oof_pred <- factor(rep(NA_character_, length(sample_ids)), levels = levels(y))
  names(oof_pred) <- sample_ids
  imp_list <- list()

  for (fi in seq_along(folds)) {
    tr <- folds[[fi]]$train_idx
    te <- folds[[fi]]$test_idx

    X_tr <- expr[, tr, drop = FALSE]
    X_te <- expr[, te, drop = FALSE]
    y_tr <- y[tr]

    # Train-only constant-gene drop (avoid constant-row issues in batch correction)
    dc <- .drop_constant_train(X_tr, X_te, eps = 1e-12)
    X_tr <- dc$X_tr; X_te <- dc$X_te

    # Fold-safe batch correction (train fit -> apply to test) if requested
    if (isTRUE(fold_batch_correction) &&
    !is.null(batch_col) && batch_col %in% names(meta)) {

  meta_tr <- meta[tr, , drop = FALSE]
  meta_te <- meta[te, , drop = FALSE]

  bc <- fold_batch_correct(
    Ytr = X_tr, Yte = X_te,
    meta_tr = meta_tr, meta_te = meta_te,
    batch_col = batch_col,
    covar_cols = if (is.null(batch_covariates)) character(0) else batch_covariates,
    label_col  = label_col
  )
  X_tr <- bc$train
  X_te <- bc$test
}

    # CV-safe preprocessing
    if (isTRUE(filter_low_expr)) {
      keep <- rowMeans(X_tr > 0) >= min_prop
      if (any(keep)) {
        X_tr <- X_tr[keep, , drop = FALSE]
        X_te <- X_te[rownames(X_tr), , drop = FALSE]
      }
    }
    if (transform == "log1p") {
      X_tr <- log1p(X_tr); X_te <- log1p(X_te)
    }
    if (isTRUE(standardize)) {
      if (!requireNamespace("matrixStats", quietly = TRUE))
        stop("matrixStats is required for standardize=TRUE.", call. = FALSE)
      mu <- matrixStats::rowMeans2(X_tr)
      sd <- matrixStats::rowSds(X_tr); sd[sd == 0] <- 1
      X_tr <- sweep(sweep(X_tr, 1, mu, "-"), 1, sd, "/")
      X_te <- sweep(sweep(X_te, 1, mu, "-"), 1, sd, "/")
    }

    # Unsupervised variance filter on TRAIN
    if (n_top > 0L && nrow(X_tr) > n_top) {
      if (!requireNamespace("matrixStats", quietly = TRUE))
        stop("matrixStats is required for variance filtering.", call. = FALSE)
      v   <- matrixStats::rowVars(X_tr)
      idx <- order(v, decreasing = TRUE)[seq_len(n_top)]
      X_tr <- X_tr[idx, , drop = FALSE]
      X_te <- X_te[rownames(X_tr), , drop = FALSE]
    }

    # Train RF
    df_tr <- data.frame(t(X_tr), check.names = FALSE)
    rf <- ranger::ranger(
      x = df_tr, y = y_tr,
      num.trees      = trees,
      importance     = importance,
      classification = TRUE,
      probability    = TRUE,
      class.weights  = class_weights,
      num.threads    = threads,
      seed           = seed + fi
    )

    # z-normalize fold importances
    imp_raw <- ranger::importance(rf)
    imp_z   <- as.numeric(scale(imp_raw))
    names(imp_z) <- names(imp_raw)
    imp_list[[fi]] <- imp_z

    # OOF predictions
    df_te <- data.frame(t(X_te), check.names = FALSE)
    pr    <- predict(rf, data = df_te)$predictions
    oof_prob[colnames(X_te), colnames(pr)] <- as.matrix(pr)
    oof_pred[colnames(X_te)] <- colnames(pr)[max.col(pr, ties.method = "first")]
  }

  # Aggregate importances
  all_genes <- unique(unlist(lapply(imp_list, names)))
  imp_mat <- do.call(cbind, lapply(imp_list, function(v) v[match(all_genes, names(v))]))
  rownames(imp_mat) <- all_genes
  imp_mat[is.na(imp_mat)] <- 0
  imp_mean <- rowMeans(imp_mat, na.rm = TRUE)
  sel_in_folds <- rowSums(imp_mat != 0)

  imp_df <- data.frame(
    gene = names(imp_mean),
    importance = as.numeric(imp_mean),
    SelectedInFolds = as.integer(sel_in_folds),
    stringsAsFactors = FALSE
  )
  imp_df <- imp_df[order(imp_df$importance, decreasing = TRUE), , drop = FALSE]

  oof <- list(prob = oof_prob, pred = oof_pred, y = y)

  # PCA variance info (for downstream PCA plots only)
  Xpca <- t(expr)  # samples x genes
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    nzv <- matrixStats::colSds(Xpca) > 0
  } else {
    nzv <- apply(Xpca, 2, sd) > 0
  }
  if (sum(nzv) >= 2) {
    pc <- stats::prcomp(Xpca[, nzv, drop = FALSE], center = TRUE, scale. = TRUE)
    var_pca <- (pc$sdev^2) / sum(pc$sdev^2)
  } else {
    var_pca <- c(0, 0)
  }

  params <- list(
    label_col = label_col,
    k = k, trees = trees, importance = importance,
    class_weights = class_weights,
    threads = threads, seed = seed,
    fold_batch_correction = fold_batch_correction,
    batch_col = batch_col,
    batch_covariates = batch_covariates,
    filter_low_expr = filter_low_expr,
    min_prop = min_prop,
    transform = transform,
    standardize = standardize,
    cv = cv,
    auto_confounds = auto_confounds
  )

  new("GeneRankFit",
      params   = params,
      oof      = oof,
      imp      = imp_df,
      features = imp_df$gene,
      final_model = NULL,
      var_pca  = var_pca,
      calibration = list())
}

