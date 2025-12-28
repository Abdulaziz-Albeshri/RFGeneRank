# ---- R/combat_helpers.R ------------------------------------------------------
# Internal helpers for ComBat with CV-safe guards.

#' Make 0-intercept model matrix; chars -> factors; drop unused levels
#' @keywords internal
#' @noRd
.to_mm0 <- function(df) {
  if (is.null(df) || ncol(df) == 0L) return(NULL)
  for (j in seq_along(df)) if (is.character(df[[j]])) df[[j]] <- factor(df[[j]])
  stats::model.matrix(~ 0 + ., data = df)
}

#' Align TEST design columns to TRAIN (add missing as zeros; order cols)
#' @keywords internal
#' @noRd
.align_mm_to_train <- function(mm_te, mm_tr) {
  if (is.null(mm_tr)) return(NULL)
  if (is.null(mm_te)) {
    mm_te <- matrix(0, nrow = 0, ncol = 0)
  }
  miss <- setdiff(colnames(mm_tr), colnames(mm_te))
  if (length(miss)) {
    mm_te <- cbind(mm_te, matrix(0, nrow = nrow(mm_te), ncol = length(miss),
                                 dimnames = list(NULL, miss)))
  }
  mm_te[, colnames(mm_tr), drop = FALSE]
}

#' Detect perfect batch~covariate confounding in TRAIN
#' - Categorical: 1-to-1 (or 1-to-*) mapping per batch row (exactly one nonzero per row
#'   and each covariate level used by at most one batch) -> confounded.
#' - Numeric: value is constant within each batch (all obs per batch identical) -> confounded.
#' @keywords internal
#' @noRd
.is_confounded_with_batch <- function(batch, covar) {
  if (is.null(covar)) return(FALSE)
  if (is.character(covar)) covar <- factor(covar)
  b <- factor(batch)

  if (is.factor(covar)) {
    tab <- table(b, covar)
    # each batch has exactly one level present
    one_per_row <- all(rowSums(tab > 0) == 1)
    # and each covar level used by at most one batch
    at_most_one_per_col <- all(colSums(tab > 0) <= 1)
    return(one_per_row && at_most_one_per_col)
  } else if (is.numeric(covar)) {
    by_batch <- tapply(covar, b, function(v) length(unique(v[!is.na(v)])))
    # every batch is constant -> effectively collinear with batch
    return(all(by_batch == 1))
  }
  FALSE
}

#' Fit ComBat on TRAIN ONLY
#' @keywords internal
#' @noRd
rfgr_combat_fit_train <- function(X_tr, batch_tr, mod_tr,
                                  par.prior = TRUE, prior.plots = FALSE) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("ComBat requires the Bioconductor package 'sva'.")
  }

  X_tr <- as.matrix(X_tr)
  b_tr <- as.vector(batch_tr)

  cb <- sva::ComBat(
    dat = X_tr,
    batch = b_tr,
    mod = mod_tr,
    par.prior = par.prior,
    prior.plots = prior.plots
  )

  mu <- rowMeans(cb, na.rm = TRUE)
  sdv <- apply(cb, 1, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1

  list(X_tr_h = cb, estimates = list(mu = mu, sd = sdv))
}

#' Apply TRAIN estimates to TEST
#' @keywords internal
#' @noRd
rfgr_combat_apply_test <- function(X_te, batch_te, mod_te, estimates) {

  X_te <- as.matrix(X_te)
  mu <- estimates$mu
  sdv <- estimates$sd

  sweep(sweep(X_te, 1, mu, "-"), 1, sdv, "/")
}


#' Fold-safe batch correction (TRAIN fit -> apply to TEST)
#' Guards:
#'  - Skip if TRAIN has <2 batches
#'  - Skip if any TEST batch is unseen in TRAIN
#'  - Drop covariates that are single-level in TRAIN
#'  - Drop covariates perfectly confounded with batch
#'  - Never include batch or label in the design
#' @keywords internal
#' @noRd
fold_batch_correct <- function(
  Ytr, Yte,
  meta_tr, meta_te,
  batch_col,
  covar_cols = character(0),
  label_col  = NULL,
  par.prior = TRUE,
  prior.plots = FALSE
) {
  b_tr <- meta_tr[[batch_col]]
  b_te <- meta_te[[batch_col]]

  # 0) TRAIN must have >=2 batches
  ub_tr <- unique(b_tr)
  if (length(ub_tr) < 2L) {
    # message("fold_batch_correct: TRAIN has <2 batches; skipping ComBat for this fold.")
    return(list(train = Ytr, test = Yte))
  }

  # 1) All TEST batches must be present in TRAIN
  ub_te <- unique(b_te)
  if (!all(ub_te %in% ub_tr)) {
    # message("fold_batch_correct: TEST has unseen batch(es); skipping ComBat for this fold.")
    return(list(train = Ytr, test = Yte))
  }

  # 2) Build safe covariate set for TRAIN
  covar_cols <- unique(setdiff(covar_cols, c(batch_col, label_col)))
  valid_covs <- character(0)
  if (length(covar_cols)) {
    keep <- vapply(covar_cols, function(cl) {
      vals <- meta_tr[[cl]]
      if (is.null(vals)) return(FALSE)
      # drop single-level factors/chars; keep numeric (even single-level numeric is dropped below if confounded)
      if (is.numeric(vals)) TRUE else (length(na.omit(unique(vals))) >= 2L)
    }, logical(1))
    valid_covs <- covar_cols[keep]
  }

  # 3) Drop covariates perfectly confounded with batch (per TRAIN)
  if (length(valid_covs)) {
    confounded <- vapply(valid_covs, function(cl) {
      .is_confounded_with_batch(b_tr, meta_tr[[cl]])
    }, logical(1))
    if (any(confounded)) valid_covs <- valid_covs[!confounded]
  }

  # 4) Build design matrices
  mod_tr <- if (length(valid_covs)) .to_mm0(meta_tr[valid_covs]) else NULL
  mod_te <- if (!is.null(mod_tr)) .align_mm_to_train(.to_mm0(meta_te[valid_covs]), mod_tr) else NULL

  # 5) Fit on TRAIN, apply to TEST
  fit <- rfgr_combat_fit_train(
    X_tr      = Ytr,
    batch_tr  = b_tr,
    mod_tr    = mod_tr,
    par.prior = par.prior,
    prior.plots = prior.plots
  )
  Yte_adj <- rfgr_combat_apply_test(
    X_te      = Yte,
    batch_te  = b_te,
    mod_te    = mod_te,
    estimates = fit$estimates
  )
  list(train = fit$X_tr_h, test = Yte_adj)
}
