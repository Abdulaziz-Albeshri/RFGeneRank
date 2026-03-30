# ---- R/id_map.R -------------------------------------------------------------

#' Map gene identifiers using AnnotationDbi
#'
#' A thin wrapper around \code{AnnotationDbi::select()} that preserves input
#' preserves order and enables straightforward conversion between identifier types (e.g., ENTREZID → SYMBOL).
#'
#' @param keys character vector of gene IDs to map (e.g., ENTREZ IDs).
#' @param from source keytype (e.g., "ENTREZID", "ENSEMBL", "SYMBOL").
#' @param to destination keytype (e.g., "SYMBOL").
#' @param OrgDb an \code{OrgDb} object. Defaults to \code{org.Hs.eg.db}.
#' @param drop_na logical; if TRUE, drop rows with missing mapped values.
#' @param unique logical; if TRUE, keep at most one mapping per input key,
#'   preferring the first match returned by \code{AnnotationDbi::select()}.
#'
#' @return data.frame with columns \code{from}, \code{to} in that order.
#' @export
#' @examples
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   x <- c("7157", "7158")  # ENTREZ IDs
#'   id_map(x, from = "ENTREZID", to = "SYMBOL")
#' }
id_map <- function(keys,
                   from,
                   to,
                   OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                   drop_na = TRUE,
                   unique = TRUE) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("AnnotationDbi is required for id_map().", call. = FALSE)
  if (length(keys) == 0) {
    return(data.frame(setNames(list(character(0)), from),
                      setNames(list(character(0)), to)))
  }

  # Run mapping
  res <- AnnotationDbi::select(
  OrgDb,
  keys   = base::unique(keys),
  keytype = from,
  columns = to
)

  # Keep only requested columns; rename consistently
  res <- res[, c(from, to), drop = FALSE]
  names(res) <- c(from, to)

  # Optionally drop NA
  if (isTRUE(drop_na)) {
    res <- res[!is.na(res[[to]]) & res[[to]] != "", , drop = FALSE]
  }

  # Optionally ensure 1:1 by taking first mapping per key
  if (isTRUE(unique) && nrow(res)) {
    first_idx <- !duplicated(res[[from]])
    res <- res[first_idx, , drop = FALSE]
  }

  # Restore original input order (for pretty joins)
  ord <- match(keys, res[[from]])
  res <- res[ord, , drop = FALSE]

  rownames(res) <- NULL
  res
}
