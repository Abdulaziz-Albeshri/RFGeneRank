#' Batch-aware cross-validation with optional ComBat modes (incl. frozen ComBat)
#'
#' @description
#' Orchestrates leakage-safe CV for RFGeneRank with two ComBat modes:
#' - "none":     no batch correction
#' - "train":    train-only batch correction with leakage-safe application to TEST
#'
#' Supports LOBO (leave-one-batch-out), Group K-Fold by batch, and standard K-Fold.
#'
#' @importFrom stats model.matrix as.formula
#'
#' @param expr numeric matrix genes x samples (continuous scale: log2CPM or log2(TPM+1))
#' @param metadata data.frame with SampleID, label_col, batch_col
#' @param label_col character; target column in metadata (e.g., "state")
#' @param batch_col character; batch/dataset column in metadata (e.g., "batch")
#' @param covariates character vector of column names to *preserve* in ComBat (added to design)
#' @param cv one of c("lobo","groupk","kfold")
#' @param k integer; number of folds for "groupk" or "kfold"
#' @param combat_mode one of c("none","train")
#' @param rf_trees integer; number of trees for ranger
#' @param seed integer; RNG seed
#' @param verbose logical; emit progress messages
#'
#' @return list(auc_by_fold, mean_auc, settings, folds_info)
#' @export
#' @examples
#'
#' # Minimal runnable example: tiny dataset, fast evaluation
#' expr <- matrix(rnorm(10 * 6), nrow = 10)
#' rownames(expr) <- paste0("gene", 1:10)
#' colnames(expr) <- paste0("sample", 1:6)
#'
#' label <- factor(rep(c("A", "B"), each = 3))
#' batch <- factor(rep(c("batch1", "batch2"), times = 3))
#'
#' metadata <- data.frame(
#'   label = label,
#'   batch = batch,
#'   row.names = colnames(expr)
#' )
#'
#' # Lightweight cross-validation: no ComBat, few trees
#' cv_res <- rfgr_crossval(
#'   expr        = expr,
#'   metadata    = metadata,
#'   label_col   = "label",
#'   batch_col   = "batch",
#'   cv          = "kfold",
#'   k           = 2,
#'   combat_mode = "none",
#'   rf_trees    = 10,
#'   verbose     = FALSE
#' )
#'
#' cv_res
rfgr_crossval <- function(expr, metadata,
                          label_col   = "state",
                          batch_col   = "batch",
                          covariates  = NULL,
                          cv          = c("lobo","groupk","kfold"),
                          k           = 5,
                          combat_mode = c("none","train"),
                          rf_trees    = 1000,
                          seed        = 1,
                          verbose     = TRUE) {

  stopifnot(is.matrix(expr), ncol(expr) == nrow(metadata))
  stopifnot(label_col %in% names(metadata))
  stopifnot(batch_col %in% names(metadata) || cv == "kfold")

  cv          <- match.arg(cv)
  combat_mode <- match.arg(combat_mode)

  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  # --------- Build folds ---------
  folds <- list(); fold_names <- character(0)

  if (cv == "lobo") {
    ub <- unique(metadata[[batch_col]])
    for (b in ub) {
      te <- which(metadata[[batch_col]] == b)
      tr <- setdiff(seq_len(ncol(expr)), te)
      folds[[length(folds) + 1L]] <- list(train_idx = tr, test_idx = te, tag = paste0("test=", b))
      fold_names <- c(fold_names, paste0("test=", b))
    }
  } else if (cv == "groupk") {
    grp <- as.factor(metadata[[batch_col]])
    ub  <- levels(grp)
    k    <- min(k, length(ub))
    # Try caret::groupKFold if available; else implement a simple splitter
    if (requireNamespace("caret", quietly = TRUE)) {
      tr_lists <- caret::groupKFold(group = grp, k = k)
      for (i in seq_along(tr_lists)) {
        tr <- tr_lists[[i]]
        te <- setdiff(seq_len(ncol(expr)), tr)
        folds[[length(folds) + 1L]] <- list(train_idx = tr, test_idx = te, tag = paste0("fold=", i))
        fold_names <- c(fold_names, paste0("fold=", i))
      }
    } else {
      .msg("Package 'caret' not installed; using a simple group split for groupK.")
      # randomize groups, then split into k roughly equal sets
      shuffled <- sample(ub)
      # chunk groups
      chunks <- split(shuffled, cut(seq_along(shuffled), breaks = k, labels = FALSE))
      for (i in seq_along(chunks)) {
        te_groups <- chunks[[i]]
        te <- which(grp %in% te_groups)
        tr <- setdiff(seq_len(ncol(expr)), te)
        folds[[length(folds) + 1L]] <- list(train_idx = tr, test_idx = te, tag = paste0("fold=", i))
        fold_names <- c(fold_names, paste0("fold=", i))
      }
    }
  } else { # kfold (sample-wise; not batch-aware)
    y  <- as.factor(metadata[[label_col]])
    k  <- min(k, length(y))
    if (requireNamespace("caret", quietly = TRUE)) {
      tr_lists <- caret::createFolds(y, k = k, returnTrain = TRUE)
      i <- 0L
      for (tr in tr_lists) {
        i <- i + 1L
        te <- setdiff(seq_len(ncol(expr)), tr)
        folds[[length(folds) + 1L]] <- list(train_idx = tr, test_idx = te, tag = paste0("fold=", i))
        fold_names <- c(fold_names, paste0("fold=", i))
      }
    } else {
      .msg("Package 'caret' not installed; using a simple unstratified K-fold split.")
      idx <- sample(seq_len(ncol(expr)))
      bins <- split(idx, cut(seq_along(idx), breaks = k, labels = FALSE))
      for (i in seq_along(bins)) {
        te <- bins[[i]]
        tr <- setdiff(seq_len(ncol(expr)), te)
        folds[[length(folds) + 1L]] <- list(train_idx = tr, test_idx = te, tag = paste0("fold=", i))
        fold_names <- c(fold_names, paste0("fold=", i))
      }
    }
  }

  # --------- Per-fold run ---------
  if (!requireNamespace("ranger", quietly = TRUE)) stop("Please install 'ranger'.")
  if (!requireNamespace("pROC", quietly = TRUE))   stop("Please install 'pROC'.")

  aucs <- numeric(length(folds))
  names(aucs) <- fold_names
  folds_info <- vector("list", length(folds)); names(folds_info) <- fold_names

  for (fi in seq_along(folds)) {
    fld <- folds[[fi]]
    tr  <- fld$train_idx; te <- fld$test_idx
    tag <- fld$tag

    X_tr <- expr[, tr, drop = FALSE]
    X_te <- expr[, te, drop = FALSE]

    meta_tr <- metadata[tr, , drop = FALSE]
    meta_te <- metadata[te, , drop = FALSE]

    # Design matrices for ComBat (preserve label + covariates)
    terms <- c(label_col, covariates)
    terms <- terms[terms %in% names(metadata)]
    form  <- if (length(terms)) reformulate(terms) else ~ 1

    mod_tr  <- stats::model.matrix(form, data = meta_tr)
    mod_te  <- stats::model.matrix(form, data = meta_te)

if (combat_mode == "train") {
  # Train-only fit + leakage-safe apply (implemented in combat_helpers.R)
  fit <- rfgr_combat_fit_train(
    X_tr     = X_tr,
    batch_tr = meta_tr[[batch_col]],
    mod_tr   = mod_tr
  )
  X_tr_h <- fit$X_tr_h
  X_te_h <- rfgr_combat_apply_test(
    X_te      = X_te,
    batch_te  = meta_te[[batch_col]],
    mod_te    = mod_te,
    estimates = fit$estimates
  )
} else {
  # "none": pass-through
  X_tr_h <- X_tr
  X_te_h <- X_te
}

    # --- Train RF and evaluate AUC ---
    y_tr <- factor(meta_tr[[label_col]])
    y_te <- factor(meta_te[[label_col]])
    if (nlevels(y_tr) != 2L || nlevels(y_te) != 2L) {
      warning(sprintf("[%s] AUC assumes binary labels; got %d/%d levels.",
                      tag, nlevels(y_tr), nlevels(y_te)))
    }

    df_tr <- data.frame(t(X_tr_h), check.names = FALSE)
    df_tr[[label_col]] <- y_tr

    fit <- ranger::ranger(
      formula      = stats::as.formula(paste(label_col, "~ .")),
      data         = df_tr,
      probability  = TRUE,
      num.trees    = rf_trees,
      respect.unordered.factors = TRUE
    )

    pr <- predict(fit, data = data.frame(t(X_te_h), check.names = FALSE))$predictions

    # Choose positive class as the second level by convention
    lv  <- levels(y_tr)
    pos <- if (length(lv) >= 2L) lv[2L] else lv[1L]

    p_pos <- if (is.matrix(pr)) {
      if (pos %in% colnames(pr)) pr[, pos] else pr[, ncol(pr)]
    } else {
      pr[, 2]
    }

    auc <- as.numeric(pROC::auc(y_te, p_pos))
    aucs[fi] <- auc
    folds_info[[fi]] <- list(tag = tag,
                             n_train = length(tr),
                             n_test  = length(te),
                             batches_train = if (batch_col %in% names(metadata)) unique(meta_tr[[batch_col]]) else NA,
                             batch_test    = if (batch_col %in% names(metadata)) unique(meta_te[[batch_col]]) else NA)
    .msg("[%s] AUC = %.3f", tag, auc)
  }

  list(
    auc_by_fold = aucs,
    mean_auc    = mean(aucs),
    settings    = list(cv = cv, k = k, combat_mode = combat_mode, rf_trees = rf_trees, seed = seed),
    folds_info  = folds_info
  )
}
