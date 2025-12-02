# =============================================================================
# align_datasets.R
# Clean, align, and merge multiple expression/metadata datasets
# =============================================================================

#' Align and merge expression + metadata (genes by intersection; strict sample match)
#'
#' @description
#' `align_datasets()` ingests parallel lists of expression tables (matrices/data.frames)
#' and metadata data.frames, cleans them, enforces strict sample alignment
#' (`rownames(metadata) == colnames(expression)`), drops metadata rows with any NA
#' (with a warning + structured report), merges by the intersection of genes, and
#' returns the merged expression, metadata, and a `SummarizedExperiment` with gene IDs
#' locked into `rowData(se)$gene_id`.
#'
#' No regex recoding is performed. If a dataset lacks a `batch` column and
#' `require_batch=TRUE`, a per-dataset batch factor ("Batch1", "Batch2", ...) is created.
#'
#' @param expr_list list of matrices or data.frames; rows = genes, cols = samples.
#'   If a data.frame contains a gene-ID column, it will be moved into rownames.
#'   All expression values are coerced to numeric (double).
#' @param meta_list list of data.frames; rows = samples (must be in rownames).
#'   Must contain a `state` column and (if `require_batch=TRUE`) a `batch` column.
#'   If `batch` is missing, it is auto-created per dataset as "Batch1", "Batch2", etc.
#' @param prefer optional selector/renamer for metadata columns. Two forms are supported:
#'   - index form: c("2:state","5:batch") (select column 2 and 5, rename to state,batch)
#'   - name map:   c("state=phenotype","batch=plate") (select phenotype, rename to state, etc.)
#'   Accepts either a single character vector applied to all datasets, or a named list specifying datasets individually.
#'   with one vector per dataset (names must match names(expr_list)/names(meta_list)).
#'   No value recoding is performed.
#' @param require_batch logical (default TRUE). If TRUE and a dataset lacks `batch`,
#'   create a per-dataset batch factor ("Batch1", "Batch2", etc.).
#' @param tag_dataset logical (default TRUE). If TRUE, add a `dataset` factor
#'   ("ds1","ds2",...) to the merged metadata.
#' @param gene_merge character(1), default "intersection". Currently only
#'   intersection is supported (keeps only genes shared by all datasets).
#' @param verbose logical. If TRUE (default), prints a brief summary and per-dataset stats.
#'
#' @return A list with:
#' \describe{
#'   \item{expr}{Merged numeric matrix (genes x samples) with syntactic, unique gene IDs.}
#'   \item{metadata}{Merged metadata (samples as rownames) with `state`/`batch` as factors.}
#'   \item{se}{A `SummarizedExperiment` with assay "expr", `colData = metadata`,
#'             and `rowData(se)$gene_id` set to current rownames.}
#'   \item{report}{List with per-dataset input/kept counts, NA-drop details, observed levels,
#'                 and final shared-gene count.}
#' }
#'
#' Align and merge expression + metadata ...
#' @examples
#' # Single toy dataset: expression matrix (genes x samples)
#' expr <- matrix(
#'   rnorm(5 * 4),
#'   nrow = 5,
#'   dimnames = list(
#'     paste0("gene", 1:5),
#'     paste0("s1_", 1:4)
#'   )
#' )
#'
#' # Matching metadata: one row per sample, with 'state' and 'batch' columns
#' meta <- data.frame(
#'   state = rep("Control", ncol(expr)),
#'   batch = rep("A",       ncol(expr))
#' )
#' rownames(meta) <- colnames(expr)
#'
#' # Lists of length 1 for expression and metadata
#' expr_list <- list(A = expr)
#' meta_list <- list(A = meta)
#'
#' aligned <- align_datasets(
#'   expr_list = expr_list,
#'   meta_list = meta_list
#' )
#' str(aligned)
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom S4Vectors DataFrame SimpleList metadata
align_datasets <- function(expr_list,
                           meta_list,
                           prefer = NULL,
                           require_batch = TRUE,
                           tag_dataset = TRUE,
                           gene_merge = "intersection",
                           verbose = TRUE) {

  # ------------------------------ validations ------------------------------
  if (!is.list(expr_list) || !is.list(meta_list))
    stop("`expr_list` and `meta_list` must be lists of equal length.")
  if (length(expr_list) != length(meta_list) || length(expr_list) < 1L)
    stop("`expr_list` and `meta_list` must be same non-zero length.")

  if (!identical(gene_merge, "intersection"))
    stop("Only gene_merge = 'intersection' is supported at the moment.")

  # dataset names
  nm <- names(expr_list)
  if (is.null(nm)) nm <- paste0("ds", seq_along(expr_list))
  if (is.null(names(meta_list))) names(meta_list) <- nm
  names(expr_list) <- nm

  # prefer handling: allow vector (apply to all) or named list
  prefer_is_list <- is.list(prefer)
  if (prefer_is_list) {
    if (is.null(names(prefer)) || !all(nm %in% names(prefer))) {
      stop("When `prefer` is a list, it must be named with dataset names: ",
           paste(nm, collapse = ", "))
    }
  }

  # ------------------------------ helpers ----------------------------------

  make_gene_ids <- function(x) {
    rn <- rownames(x)
    if (is.null(rn) || any(!nzchar(rn))) {
      # try first column as IDs if character-like and mostly unique
      if (is.data.frame(x)) {
        char_cols <- vapply(x, function(col) !is.numeric(col), logical(1))
        if (any(char_cols)) {
          j <- which(char_cols)[1]
          ids <- as.character(x[[j]])
          if (!anyNA(ids) && length(unique(ids)) / length(ids) > 0.99) {
            rn <- ids
            x <- x[, setdiff(colnames(x), colnames(x)[j]), drop = FALSE]
          }
        }
      }
    }
    if (is.null(rn) || any(!nzchar(rn))) {
      stop("Gene rownames are missing and no suitable ID column was detected.")
    }
     rn <- as.character(rn)
    rn <- trimws(rn)

    # Fix R's "X" prefix on numeric IDs (X100287102 -> 100287102)
    rn <- sub("^X(?=\\d+$)", "", rn, perl = TRUE)

    # Ensure uniqueness without mangling the base IDs
    rn <- make.unique(rn, sep = "_dup")

    list(object = x, rownames = rn)
  }

  as_numeric_matrix <- function(x, rn) {
    if (!is.matrix(x)) x <- as.matrix(x)
    Xn <- matrix(as.numeric(x), nrow = nrow(x),
             dimnames = list(rn, colnames(x)))
    if (any(!is.finite(Xn)))
      stop("Non-finite values detected in expression matrix after coercion.")
    storage.mode(Xn) <- "double"
    Xn
  }

  apply_prefer <- function(M, pref, ds_name) {
    if (is.null(pref) || length(pref) == 0) return(M)
    out <- list(); newnames <- character(0)

    # support "idx:name" and "new=old"
    parse_one <- function(token) {
      if (grepl(":", token, fixed = TRUE)) {
        p <- strsplit(token, ":", fixed = TRUE)[[1]]
        if (length(p) != 2) stop("Invalid prefer token '", token,
                                 "'. Use '2:state' or 'state=phenotype'.")
        idx <- as.integer(p[1])
        if (is.na(idx) || idx < 1 || idx > ncol(M))
          stop("Invalid column index in prefer for dataset ", ds_name, ": ", token)
        list(name = p[2], vec = M[[idx]])
      } else if (grepl("=", token, fixed = TRUE)) {
        p <- strsplit(token, "=", fixed = TRUE)[[1]]
        if (length(p) != 2) stop("Invalid prefer token '", token,
                                 "'. Use '2:state' or 'state=phenotype'.")
        if (!p[2] %in% colnames(M))
          stop("Column '", p[2], "' not found in metadata for dataset ", ds_name, ".")
        list(name = p[1], vec = M[[p[2]]])
      } else {
        stop("Unrecognized prefer token '", token, "'.")
      }
    }

    for (tok in pref) {
      pr <- parse_one(tok)
      out[[length(out) + 1L]] <- pr$vec
      newnames <- c(newnames, pr$name)
    }
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    colnames(out) <- newnames
    out
  }

  # warning utility (also returns a structured list for the report)
  warn_na_drops <- function(dataset_name, dropped_ids, affected_cols) {
    N <- length(dropped_ids)
    if (N == 0) return(list(n_dropped = 0L, dropped_ids = character(0), affected_cols = character(0)))
    msg <- paste0(
      "Dropped ", N, " sample(s) from dataset ", dataset_name,
      " due to missing metadata values.\n",
      "  * Affected columns included: ", paste(unique(affected_cols), collapse = ", "), "\n",
      "  * Example dropped sample IDs: ",
      paste(head(dropped_ids, 5), collapse = ", "),
      if (N > 5) ", ..." else ""
    )
    warning(msg, call. = FALSE)
    list(n_dropped = N, dropped_ids = dropped_ids, affected_cols = unique(affected_cols))
  }

  # --------------------------- per-dataset pass -----------------------------
  per_ds <- setNames(vector("list", length(nm)), nm)
  report  <- list(per_dataset = setNames(vector("list", length(nm)), nm))

  for (i in seq_along(nm)) {
    ds <- nm[i]
    X0 <- expr_list[[i]]
    M0 <- as.data.frame(meta_list[[i]], stringsAsFactors = FALSE)

    if (nrow(M0) == 0) stop("Metadata '", ds, "' is empty.")
    if (is.null(rownames(M0)) || any(!nzchar(rownames(M0))))
      stop("Metadata '", ds, "' must have rownames set to sample IDs (no `sample_col` fallback).")

    # ---- Expression: gene IDs + numeric coercion
    g <- make_gene_ids(X0)
    X1 <- g$object
    rn <- g$rownames
    X  <- as_numeric_matrix(X1, rn)

    # ---- Metadata: prefer (select/rename) without recoding
    pref_i <- if (prefer_is_list) prefer[[ds]] else prefer
    M <- if (!is.null(pref_i)) apply_prefer(M0, pref_i, ds) else M0

    # ---- Require state; handle/require batch
    if (!("state" %in% colnames(M)))
      stop("Dataset '", ds, "': metadata must contain 'state' (use `prefer` to rename if needed).")
    if (require_batch && !("batch" %in% colnames(M))) {
      M$batch <- paste0("Batch", i)  # one level for this dataset
    }

    # ---- Strict NA policy: drop any row with any NA (warn + report)
    na_mask <- !stats::complete.cases(M)
    drop_ids <- rownames(M)[na_mask]
    aff_cols <- names(M)[colSums(is.na(M)) > 0]
    drop_info <- warn_na_drops(ds, drop_ids, aff_cols)
    if (any(na_mask)) M <- M[!na_mask, , drop = FALSE]

    # ---- Factors for state/batch (batch may now exist)
    M$state <- factor(M$state)
    if ("batch" %in% colnames(M)) M$batch <- factor(M$batch)

    # ---- Strict alignment: exact IDs only
    ids_expr <- colnames(X)
    ids_meta <- rownames(M)
    overlap  <- intersect(ids_expr, ids_meta)
    if (length(overlap) == 0L)
      stop("Dataset '", ds, "': no common sample IDs between expression columns and metadata rownames.")

    # enforce exact match/order
    X <- X[, overlap, drop = FALSE]
    M <- M[ overlap, , drop = FALSE]
    stopifnot(identical(colnames(X), rownames(M)))

    # ---- per-dataset report
    report$per_dataset[[ds]] <- list(
      genes_in   = nrow(X),
      samples_in = length(ids_expr),
      samples_kept = ncol(X),
      na_drop = drop_info,
      levels = list(
        state = levels(M$state),
        batch = if ("batch" %in% colnames(M)) levels(M$batch) else character(0)
      )
    )

    per_ds[[ds]] <- list(expr = X, meta = M)
  }

    # ------------------------------ merge pass -------------------------------
   # ------------------------------ merge pass -------------------------------
  # gene intersection (safe; no NA)
  shared <- Reduce(intersect, lapply(per_ds, function(z) rownames(z$expr)))
  if (length(shared) == 0L)
    stop("No shared genes across datasets after cleaning.")

  n_ds <- length(per_ds)

  if (n_ds == 1L) {
    # --- Single-dataset fast path (robust) ---
    X_merged <- per_ds[[1]]$expr[shared, , drop = FALSE]
    M_merged <- per_ds[[1]]$meta

    if (isTRUE(tag_dataset)) {
      M_merged$dataset <- factor("ds1", levels = "ds1")
    }

  } else {
    # --- Multi-dataset path ---
    X_merged <- do.call(cbind, lapply(per_ds, function(z) z$expr[shared, , drop = FALSE]))
    M_merged <- do.call(rbind,  lapply(per_ds, function(z) z$meta))

    if (isTRUE(tag_dataset)) {
      ds_tags <- rep(names(per_ds),
                     times = vapply(per_ds, function(z) ncol(z$expr), integer(1)))
      M_merged$dataset <- factor(ds_tags, levels = unique(ds_tags))
    }
  }

  # Force final metadata row order to match expression column order
  M_merged <- M_merged[colnames(X_merged), , drop = FALSE]

  # final invariants (with diagnostics)
  if (!identical(colnames(X_merged), rownames(M_merged))) {
    bad <- which(colnames(X_merged) != rownames(M_merged))
    msg <- paste0(
      "Final alignment failed: colnames(expr) != rownames(metadata).\n",
      "  length(expr cols) = ", ncol(X_merged), ", length(meta rows) = ", nrow(M_merged), "\n",
      "  First mismatches at positions: ", paste(head(bad, 10), collapse = ", "), "\n",
      "  Example expr IDs: ", paste(head(colnames(X_merged)[bad], 5), collapse = ", "), "\n",
      "  Example meta IDs: ", paste(head(rownames(M_merged)[bad], 5), collapse = ", ")
    )
    stop(msg)
  }

  # ensure state/batch are factors in merged meta
  if ("state" %in% names(M_merged)) M_merged$state <- factor(M_merged$state)
  if ("batch" %in% names(M_merged)) M_merged$batch <- factor(M_merged$batch)

  # ------------------------- SummarizedExperiment --------------------------
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Please install 'SummarizedExperiment' and 'S4Vectors'.")
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = S4Vectors::SimpleList(expr = X_merged),
    colData = S4Vectors::DataFrame(M_merged)
  )
  # lock gene IDs
  SummarizedExperiment::rowData(se)$gene_id <- rownames(se)

  # ------------------------------- summary ---------------------------------
  if (isTRUE(verbose)) {
    message(sprintf("[align_datasets] merged %d dataset(s)\n", length(per_ds)))
    message(sprintf("Shared genes: %d | Samples: %d\n", nrow(X_merged), ncol(X_merged)))
  }

  report$global <- list(
    datasets      = names(per_ds),
    shared_genes  = length(shared),
    samples_final = ncol(X_merged)
  )

  list(
    expr     = X_merged,
    metadata = M_merged,
    se       = se,
    report   = report
  )
}
