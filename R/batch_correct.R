# ---- R/batch_correct.R ------------------------------------------------------
# Robust, fold-safe batch correction.

# Build design with intercept + optional covariates + batch (treatment coding).
# Automatically drops empty levels in the train fold.
.build_design <- function(meta, batch_col, covar_cols = NULL) {
  if (is.null(batch_col) || !batch_col %in% names(meta))
    stop("batch_col not found in metadata.", call. = FALSE)

  df <- meta
  df[[batch_col]] <- factor(df[[batch_col]])  # ensure factor
  # keep only requested covariates that exist
  keep_cov <- intersect(covar_cols %||% character(0), names(df))
  # Build formula: ~ 1 + covars + batch
  if (length(keep_cov)) {
    fml <- stats::as.formula(
      paste("~ 1 +", paste(keep_cov, collapse = " + "), "+", batch_col)
    )
  } else {
    fml <- stats::as.formula(paste("~ 1 +", batch_col))
  }
  X <- stats::model.matrix(fml, data = df)

  # Identify which columns correspond to batch indicators.
  # model.matrix() with treatment coding will produce columns that start with the batch name.
  batch_cols <- grep(paste0("^", batch_col), colnames(X))
  attr(X, "batch_cols") <- batch_cols
  X
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Robust linear solve: (X'X + lambda I)^{-1} X'Y with fallback to pseudoinverse
.solve_beta <- function(Xtr, Ytr, lambda = 1e-6) {
  Xt  <- t(Xtr)
  XtX <- Xt %*% Xtr
  # ridge
  p <- ncol(Xtr)
  XtX_reg <- XtX + diag(lambda, p)
  XtY <- Xt %*% t(Ytr)     # p x G

  beta <- tryCatch(solve(XtX_reg, XtY),
                   error = function(e) {
                     # Final fallback: pseudoinverse
                     if (!requireNamespace("MASS", quietly = TRUE))
                       stop("MASS is required for pseudoinverse fallback.", call. = FALSE)
                     MASS::ginv(XtX) %*% XtY
                   })
  beta  # p x G
}

# Subtract ONLY batch contribution from Y (genes x samples)
.lm_apply_remove_batch <- function(Y, X, beta, batch_cols) {
  if (length(batch_cols) == 0) return(Y)  # nothing to remove
  keep <- matrix(0, nrow = nrow(beta), ncol = ncol(beta),
                 dimnames = dimnames(beta))
  keep[batch_cols, ] <- beta[batch_cols, , drop = FALSE]
  batch_pred <- X %*% keep            # samples x genes
  Y - t(batch_pred)                   # genes x samples
}

# Public: fold-safe correction on train+test matrices
fold_batch_correct <- function(Ytr, Yte, meta_tr, meta_te,
                               batch_col, covar_cols = NULL) {
  # Build train design; drops empty batch levels automatically
  Xtr <- .build_design(meta_tr, batch_col, covar_cols)

  # If train has no batch columns (single batch level after treatment coding), skip correction
  bcols_tr <- attr(Xtr, "batch_cols")
  if (length(bcols_tr) == 0) {
    return(list(train = Ytr, test = Yte))
  }

  # Build test design using same formula; align its columns to train design
  Xte_full <- .build_design(meta_te, batch_col, covar_cols)
  # Add any missing columns present in train design
  miss_cols <- setdiff(colnames(Xtr), colnames(Xte_full))
  if (length(miss_cols)) {
    add <- matrix(0, nrow = nrow(Xte_full), ncol = length(miss_cols),
                  dimnames = list(rownames(Xte_full), miss_cols))
    Xte_full <- cbind(Xte_full, add)
  }
  # Reorder exactly like train
  Xte <- Xte_full[, colnames(Xtr), drop = FALSE]

  # Fit beta on TRAIN only (ridge-regularized; pseudoinverse fallback)
  beta <- .solve_beta(Xtr, Ytr, lambda = 1e-6)

  # Apply: subtract batch contribution with train-learned coefficients
  Ytr_corr <- .lm_apply_remove_batch(Ytr, Xtr, beta, bcols_tr)
  Yte_corr <- .lm_apply_remove_batch(Yte, Xte, beta, bcols_tr)

  list(train = Ytr_corr, test = Yte_corr)
}
