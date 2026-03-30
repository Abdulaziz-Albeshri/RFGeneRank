# ---- R/top_genes.R ----------------------------------------------------------

#' Extract top predictive genes from a GeneRankFit
#'
#' Returns a ranked list from the aggregated, fold-normalized importance
#' stored in the fitted object. Optionally adds an ID mapping column (e.g.,
#' ENTREZID -> SYMBOL) using \code{\link{id_map}}.
#'
#' @param fit a \code{GeneRankFit} object (from \code{rank_genes()}).
#' @param n integer; number of top genes to return (default 100).
#' @param map logical; if TRUE, map \code{from -> to} and add a \code{mapped} column.
#' @param OrgDb an \code{OrgDb} object for mapping (default \code{org.Hs.eg.db}).
#' @param from source keytype used for the gene identifiers (e.g., "ENTREZID", "SYMBOL", "ENSEMBL").
#' @param to destination keytype (e.g., "SYMBOL").
#'
#' @return 
#' A list with:
#' \describe{
#'   \item{gene}{character vector of top gene IDs in \code{from} keytype}
#'   \item{table}{data.frame with columns \code{gene}, \code{importance},
#'                \code{SelectedInFolds}, and optional \code{mapped}}
#' }
#'
#' @export
#' @examples
#' gene_scores <- data.frame(
#'   gene       = paste0("gene", 1:6),
#'   importance = c(5, 4, 3, 2, 1, 0)
#' )
#'
#' fit <- methods::new("GeneRankFit", imp = gene_scores)
#'
#' tg <- top_genes(fit, n = 3)
#' tg$gene
#' tg$table
#' # In practice, top_genes() is used on a GeneRankFit object.
#' # For example, after running a full RFGeneRank pipeline:
#' #
#' #   fit <- gene_rank(se, genes = rownames(se), ...)
#' #   head(top_genes(fit, n = 20))
#' #
#' # where 'fit' is a GeneRankFit containing gene importance scores.
top_genes <- function(fit,
                      n = 100,
                      map = FALSE,
                      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                      from = "SYMBOL",
                      to   = "SYMBOL") {
  if (!inherits(fit, "GeneRankFit"))
    stop("fit must be a GeneRankFit.", call. = FALSE)

  imp <- .rfgr_imp(fit)
  if (!is.data.frame(imp) || !all(c("gene", "importance") %in% names(imp))) {
    stop("Stored importance table is missing required columns 'gene' and 'importance'.", call. = FALSE)
  }

  # Order by importance (already aggregated across folds) and take top n
  df <- imp[order(imp$importance, decreasing = TRUE), , drop = FALSE]
  if (nrow(df) == 0L) {
    return(list(gene = character(0), table = df))
  }
  df <- head(df, n) # return only the top n ranked genes

  # Optional ID mapping
  if (isTRUE(map)) {
    mp <- id_map(df$gene, from = from, to = to, OrgDb = OrgDb,
                 drop_na = FALSE, unique = TRUE)
    # build a mapped vector preserving order; if any missing, keep NA
    mapped_vec <- setNames(mp[[to]], mp[[from]])[df$gene]
    df$mapped <- unname(mapped_vec)
  }

  # Return both the vector and a richer table
  list(
    gene  = df$gene,
    table = df
  )
}
