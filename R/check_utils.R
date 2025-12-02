# ---- R/check_utils.R --------------------------------------------------------

.check_is_matrix <- function(m, name) {
  if (!is.matrix(m)) stop(name, " must be a matrix [genes x samples].", call. = FALSE)
  if (!is.numeric(m)) stop(name, " must be numeric (counts/TPM/CPM).", call. = FALSE)
  if (nrow(m) < 2L || ncol(m) < 2L) stop(name, " must have >=2 genes and >=2 samples.", call. = FALSE)
}

.check_meta <- function(meta, label_col, batch_col = NULL) {
  if (!is.data.frame(meta)) stop("metadata must be a data.frame.", call. = FALSE)
  if (!label_col %in% names(meta)) stop("metadata is missing label_col: ", label_col, call. = FALSE)
  if (!is.null(batch_col) && !batch_col %in% names(meta)) stop("metadata is missing batch_col: ", batch_col, call. = FALSE)
}

.check_same_samples <- function(expr, meta) {
  common <- intersect(colnames(expr), rownames(meta))
  if (length(common) < 4L) stop("Fewer than 4 overlapping samples between expression and metadata.", call. = FALSE)
  expr <- expr[, common, drop = FALSE]
  meta <- meta[common, , drop = FALSE]
  list(expr = expr, meta = meta)
}

.message_once <- function(...) message(...)

# Lightweight seed helper (works even if users set RNG elsewhere)
.set_seed <- function(seed) {
  invisible(NULL)
}
