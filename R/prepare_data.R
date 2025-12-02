# ---- R/prepare_data.R -------------------------------------------------------
#' Prepare matrices + metadata: align, (log) transform, batch correct, prefilter
#' @param mats list of numeric matrices (genes-by-samples)
#' @param metas list of data.frames (rownames = sample IDs)
#' @param label_col outcome column in metadata
#' @param batch_col batch column in metadata (or NULL)
#' @param n_var keep top-N variable genes globally (unsupervised)
#' @param log1p logical; log1p transform before variance filter
#' @param batch_method "none","combat","limma","combat_seq"
#' @param counts logical; TRUE if raw counts (for combat_seq)
#' @param filter_in_cv logical; if TRUE, skip global variance filter and let CV do it inside folds
#' @param batch_correction_scope "global" or "fold" (fold correction happens inside CV)
#' @param batch_covariates optional character vector of metadata column names used as covariates in batch correction
#' @return SummarizedExperiment
#' @export
#' @examples
#' set.seed(1)
#'
#' expr <- matrix(stats::rnorm(20 * 10), nrow = 20)
#' rownames(expr) <- paste0("gene", 1:20)
#' colnames(expr) <- paste0("sample", 1:10)
#'
#' label  <- factor(rep(c("A", "B"), each = 5))
#' batch  <- factor(rep(c("batch1", "batch2"), times = 5))
#'
#' # Build lists of matrices + metadata as expected by prepare_data()
#' mats  <- list(expr)
#' metas <- list(data.frame(
#'   label = label,
#'   batch = batch,
#'   row.names = colnames(expr)
#' ))
#'
#' prep <- prepare_data(
#'   mats      = mats,
#'   metas     = metas,
#'   label_col = "label",
#'   batch_col = "batch"
#' )
#' prep
prepare_data <- function(
  mats, metas, label_col, batch_col = NULL,
  n_var = 5000, log1p = TRUE,
  batch_method = c("none","combat","limma","combat_seq"),
  counts = FALSE,
  filter_in_cv = FALSE,
  batch_correction_scope = c("global","fold"),
  batch_covariates = NULL
) {
  batch_method <- match.arg(batch_method)
  batch_correction_scope <- match.arg(batch_correction_scope)

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("SummarizedExperiment is required.", call. = FALSE)

  expr <- mats[[1]]; meta <- metas[[1]]
  .check_is_matrix(expr, "expr"); .check_meta(meta, label_col, batch_col)

  if (is.null(rownames(meta))) {
    id_col <- if ("sample" %in% names(meta)) "sample" else names(meta)[1]
    rownames(meta) <- meta[[id_col]]
  }
  aligned <- .check_same_samples(expr, meta)
  expr <- aligned$expr; meta <- aligned$meta

  if (isTRUE(log1p)) expr <- log1p(expr)

  if (!isTRUE(filter_in_cv) && n_var > 0 && nrow(expr) > n_var) {
    v <- matrixStats::rowVars(expr)
    idx <- order(v, decreasing = TRUE)[seq_len(n_var)]
    expr <- expr[idx, , drop = FALSE]
  }

  # Global batch correction is applied here (train/test leakage is acceptable by design); fold-safe correction is instead performed within cross-validation.
  if (batch_method != "none" && !is.null(batch_col) && batch_correction_scope == "global") {
    b <- factor(meta[[batch_col]])
    if (batch_method == "combat") {
      if (!requireNamespace("sva", quietly = TRUE)) stop("sva required for ComBat.", call. = FALSE)
      mod <- NULL  # no outcome in design
      expr <- sva::ComBat(dat = as.matrix(expr), batch = b, mod = mod,
                          par.prior = TRUE, prior.plots = FALSE)
    } else if (batch_method == "limma") {
      if (!requireNamespace("limma", quietly = TRUE)) stop("limma required.", call. = FALSE)
      # optional covariates (but no outcome)
      design <- if (!is.null(batch_covariates) && length(batch_covariates)) {
        stats::model.matrix(~ ., data = meta[, intersect(batch_covariates, names(meta)), drop = FALSE])
      } else stats::model.matrix(~ 1, data = meta)
      expr <- limma::removeBatchEffect(expr, batch = b, design = design)
    } else if (batch_method == "combat_seq") {
  if (!requireNamespace("sva", quietly = TRUE)) stop("sva required for ComBat-Seq.", call. = FALSE)

  # ComBat-Seq expects *integer counts*. Build a safe 'cnt' matrix:
  # - If data were log1p(counts), invert then round/clamp;
  # - If not logged, round/clamp directly;
  # - If 'counts' is FALSE, still coerce (warned).
  if (!counts) .message_once("ComBat-Seq requested but counts=FALSE; coercing to integers.")
  cnt <- if (isTRUE(log1p)) exp(expr) - 1 else expr
  cnt <- pmax(round(cnt), 0L)  # integers, no negatives

  # Run ComBat-Seq on integer counts
  expr_adj <- sva::ComBat_seq(as.matrix(cnt), batch = b, covar_mod = NULL, full_mod = FALSE)

  # Re-apply log1p if it was previously applied in upstream processing
  expr <- if (isTRUE(log1p)) log1p(expr_adj) else expr_adj


  }

}
    # keep covariates info in colData so rank_genes() can do fold-safe correction if requested
  meta$..batch_col <- if (!is.null(batch_col)) batch_col else NA_character_
  meta$..batch_covariates <- I(list(batch_covariates))

     # ---- Confounding diagnostics (messages only; no behavior change) ----
  if (!is.null(batch_col) && batch_col %in% colnames(meta) &&
      !is.null(label_col) && label_col %in% colnames(meta)) {

    lab <- meta[[label_col]]
    bat <- meta[[batch_col]]

    if (length(unique(lab)) > 1 && length(unique(bat)) > 1) {
      tab <- table(lab, bat)

      # Empty cells flag (hard confounding)
      if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
        warning("Empty cells in labelxbatch; potential confounding.")
      }

      # chi-square association (strong association warning)
      pchi <- try(stats::chisq.test(tab), silent = TRUE)
      if (!inherits(pchi, "try-error") && is.finite(pchi$p.value) && pchi$p.value < 1e-4) {
        warning("Strong label-batch association (chi-square p<1e-4). ",
                "Consider protecting label during batch correction (train-only design).")
      }

      # Small batches (unstable estimates)
      bcounts <- colSums(tab)
      if (any(bcounts < 5)) {
        warning("Some batches have <5 samples; batch estimates may be unstable.")
      }
    }
  }

  # ---- Build SE and annotate metadata so rank_genes() can enforce "no double correction" ----
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(expr = as.matrix(expr)),
    colData = meta
  )

  # metadata getter/setter are in S4Vectors
  md <- S4Vectors::metadata(se)
  if (is.null(md)) md <- list()

  # Keep defaults for downstream (optional)
  md$`..batch_col`        <- batch_col
  md$`..batch_covariates` <- batch_covariates

  # Mark whether a GLOBAL (one-shot) correction was actually applied here
  if (!is.null(batch_col) &&
      identical(batch_correction_scope, "global") &&
      !identical(batch_method, "none")) {
    md$`..batch_corrected`        <- TRUE
    md$`..batch_correction_scope` <- "global"
  } else {
    md$`..batch_corrected`        <- FALSE
    md$`..batch_correction_scope` <- if (identical(batch_correction_scope, "fold")) "fold" else "none"
  }

  S4Vectors::metadata(se) <- md  # setter (must assign back)
  return(se)
}
