# ---- R/confounding.R ---------------------------------------------------------
#' Detect batch ~ label association and propose next steps (internal)
#'
#' Runs a chi-squared test on the batch-by-label table, reports Cramer's V, and
#' returns actionable recommendations for prepare_data(), rank_genes(), and
#' rfgr_crossval() depending on severity and dataset structure.
#'
#' @keywords internal
#' @noRd
.check_confounding <- function(labels, batch, alpha = 0.01, v_thresh = 0.20) {
  out <- list(
    # statistics
    p = NA_real_, cramerV = NA_real_, warn = FALSE, severity = "none", msg = NULL,
    # context
    n_batches = NA_integer_, class_counts = NULL, per_batch_props = NULL,
    single_class_batches = character(0),
    # machine-readable suggestions for the pipeline
    recommendations = list(
      prepare_data  = NULL,   # list(require_batch, tag_dataset, batch_method)
      rank_genes    = NULL,   # list(cv, k, fold_batch_correction, batch_col, n_top, trees, importance, use_class_weights)
      rfgr_crossval = NULL,   # list(cv, combat_mode, rf_trees, seed)
      notes         = character(0)
    )
  )

  # Basic guards
  if (is.null(labels) || is.null(batch) || anyNA(labels) || anyNA(batch)) return(out)
  labels <- as.factor(labels); batch <- as.factor(batch)
  tab <- table(batch, labels)
  if (any(dim(tab) < 2)) return(out)

  # Chi-squared & Cramer's V
  cs <- tryCatch(
  stats::chisq.test(tab),
  warning = function(w) stats::chisq.test(tab, simulate.p.value = TRUE),
  error   = function(e) NA
)
  n  <- sum(tab)
  k  <- min(nrow(tab), ncol(tab))
  V  <- sqrt((as.numeric(cs$statistic) / n) / (k - 1))

  # Severity buckets (thresholds can be adjusted based on analytical requirements)
  sev <- if (is.nan(V)) "none" else if (cs$p.value < alpha && V >= 0.35) {
    "severe"
  } else if (cs$p.value < alpha && V >= v_thresh) {
    "moderate"
  } else if (cs$p.value < 0.05 && V >= 0.10) {
    "mild"
  } else "none"

  # Per-batch proportions and single-class batches
  props <- prop.table(tab, 1)
  single_class <- rownames(tab)[apply(props, 1, function(x) any(x == 1))]
  class_counts <- colSums(tab)

  # Populate core outputs
  out$p <- cs$p.value
  out$cramerV <- V
  out$warn <- sev %in% c("moderate","severe")
  out$severity <- sev
  out$n_batches <- nrow(tab)
  out$class_counts <- class_counts
  out$per_batch_props <- props
  out$single_class_batches <- single_class
  out$msg <- sprintf(
    "Batch~label: p=%.3g, Cramer's V=%.3f%s",
    out$p, out$cramerV,
    if (out$warn) paste0(" (", sev, " association)") else ""
  )

  # ---------------------- Recommendations -------------------------------------
  # prepare_data(): avoid global ComBat when confounded; always track batches
  rec_prepare <- list(
    require_batch = TRUE,
    tag_dataset   = TRUE,
    batch_method  = if (out$warn) "none" else "none" # keep none before CV; correction happens in-fold
  )

  # Choose CV mode
  cv_mode <- if (out$warn) "lobo" else "kfold"
  if (!out$warn && out$n_batches >= 4) cv_mode <- "groupk"
  if (length(single_class)) cv_mode <- "lobo"  # strict if any batch is single-class

  # Class imbalance heuristic (>=1.5x)
  imb_ratio <- {
    cc <- as.numeric(class_counts)
    if (length(cc) >= 2) max(cc) / min(cc) else 1
  }
  use_class_weights <- imb_ratio >= 1.5

  # Cross-val ComBat mode: frozen when confounded / single-class batches
  combat_mode <- if (out$warn || length(single_class)) "train" else "none"

  rec_rank <- list(
    cv                      = cv_mode,          # if rank_genes() supports cv=
    k                       = 5L,               # ignored by LOBO
    fold_batch_correction   = TRUE,             # train-only correction in folds
    batch_col               = "batch",
    n_top                   = 500L,
    trees                   = 1000L,
    importance              = "permutation",
    use_class_weights       = use_class_weights
  )

  rec_cv <- list(
    cv          = cv_mode,       # "lobo" | "groupk" | "kfold"
    combat_mode = combat_mode,   # "train" (frozen) when confounded
    rf_trees    = 1000L,
    seed        = 2025L
  )

  notes <- character(0)
  if (length(single_class)) {
    notes <- c(notes, sprintf(
      "Single-class batches: %s (prefer LOBO; consider using these as external test only).",
      paste(single_class, collapse = ", ")
    ))
  }
  if (use_class_weights) {
    notes <- c(notes, "Class imbalance detected; enable class weights.")
  }

  out$recommendations <- list(
    prepare_data  = rec_prepare,
    rank_genes    = rec_rank,
    rfgr_crossval = rec_cv,
    notes         = notes
  )

  out
}
