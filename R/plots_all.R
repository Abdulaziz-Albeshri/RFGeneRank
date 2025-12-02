# ---- R/plots_all.R -----------------------------------------------------------
# Unified plotting API for RFGeneRank.
# All public plotting functions live here to simplify maintenance & docs.

# small helper (internal)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- helpers used by plot_shap_dependence ------------------------------------
.make_age_numeric <- function(cd, age_col) { ... }
.age_breaks <- function(x, width = 10, start = NULL) { ... }
.parse_bin  <- function(lbl) { ... }
.summarize_age_bins <- function(df, width = 10, alpha = 0.05) { ... }
.add_band_highlight <- function(p, bins_tbl, alpha = 0.12,
                                mode = c("significant","top_k"), top_k = 1) { ... }

.rfgr_sex_colors <- c(Female="#ff5ca8", Male="#00bcd4",
                      female="#ff5ca8", male="#00bcd4",
                      unknown="#7aa6ff", all="grey40")

# ==============================================================================
# 1) Core plots derived from the primary working implementation
# ==============================================================================

#' Feature importance (top genes)
#'
#' Uses fold-normalized mean importances aggregated across CV folds.
#' To display gene SYMBOLs on the axes, set \code{map_to_symbol = TRUE}.
#' (requires org.Hs.eg.db and RFGeneRank::top_genes()).
#'
#' @param fit GeneRankFit from \code{rank_genes()}.
#' @param top Integer; number of genes to display (default 30).
#' @param map_to_symbol Logical; map gene IDs to symbols if available (default FALSE).
#' @param from Keytype for input IDs (default "ENTREZID").
#' @param to   Keytype for output labels (default "SYMBOL").
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point labs theme_classic coord_flip scale_x_discrete
#' @examples
#' # Toy expression matrix: genes x samples
#' expr <- matrix(
#'   rnorm(10 * 20),
#'   nrow = 10,
#'   dimnames = list(
#'     paste0("gene", 1:10),
#'     paste0("sample", 1:20)
#'   )
#' )
#'
#' # Binary phenotype stored in 'label' column
#' y <- factor(rep(c("Control", "Case"), each = 10))
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(expr = expr),
#'   colData = data.frame(
#'     label    = y,
#'     row.names = colnames(expr)
#'   )
#' )
#'
#' # Fit a small GeneRank model (note: no 'genes' argument)
#' fit <- rank_genes(
#'   se        = se,
#'   label_col = "label",
#'   n_top     = 10,
#'   trees     = 200
#' )
#'
#' # Plot feature importance for the top-ranked genes
#' plot_importance(fit)
plot_importance <- function(fit, top = 30,
                            map_to_symbol = FALSE,
                            from = "ENTREZID", to = "SYMBOL") {
  stopifnot(inherits(fit, "GeneRankFit"))
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plotting.", call. = FALSE)

  # top rows from importance table
  df <- fit@imp[order(fit@imp$importance, decreasing = TRUE), , drop = FALSE]
  df <- head(df, top)
  df$gene <- factor(df$gene, levels = rev(df$gene))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene, y = importance)) +
    ggplot2::geom_point() +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = sprintf("Top %d predictive genes", top),
      x = "Gene", y = "Fold-normalized importance"
    )

  if (isTRUE(map_to_symbol)) {
    # Use top_genes() to map IDs -> SYMBOLs for labels only
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      warning("org.Hs.eg.db is not installed; showing original IDs.", call. = FALSE)
      return(p)
    }
    if (!"top_genes" %in% getNamespaceExports("RFGeneRank")) {
      warning("top_genes() not available; showing original IDs.", call. = FALSE)
      return(p)
    }
    tg <- RFGeneRank::top_genes(fit, n = top, map = TRUE, from = from, to = to)$table
    labs <- setNames(tg$mapped, tg$gene) # ensure mapping follows plot order
    p <- p + ggplot2::scale_x_discrete(labels = labs)
  }

  p
}

#' Signed feature importance (directional effect)
#'
#' Visualizes direction-aware importances produced by \code{sign_importance()}.
#' If \code{tab} is NULL, the function will try to read \code{fit@imp} and
#' require a \code{signed_importance} column to be present there.
#'
#' @param fit GeneRankFit (optional if \code{tab} is supplied)
#' @param tab data.frame from \code{sign_importance()} with columns:
#'            gene, importance, direction (-1/0/1), signed_importance
#' @param top Integer; number of genes to show (default 30)
#' @param map_to_symbol Logical; map x-axis labels to SYMBOLs (default FALSE)
#' @param from Keytype for input IDs (default "ENTREZID")
#' @param to   Keytype for output labels (default "SYMBOL")
#' @param show_legend Logical; show legend (default TRUE)
#' @param palette Optional named vector for fill colors, e.g.
#'   \code{c(`-1`="#3182bd", `1`="#de2d26", `0`="#9e9e9e")}
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_col labs theme_classic coord_flip scale_x_discrete scale_fill_manual
#' @examples
#' # Toy signed importance table for 5 genes
#' signed_imp <- data.frame(
#'   gene        = paste0("gene", 1:5),
#'   importance  = c(0.5, 0.4, 0.3, 0.2, 0.1),
#'   signed_effect = c(0.5, -0.4, 0.3, -0.2, 0.1)
#' )
#'
#' head(signed_imp)
#'
#' \donttest{
#' # In practice, plot_sign_importance() is used on a GeneRankFit object
#' # after computing signed importances, for example:
#' #
#' #   fit <- gene_rank(se, genes = rownames(se), ...)
#' #   sig_imp <- sign_importance(fit, X = expr_matrix, y = outcome)
#' #   plot_sign_importance(fit, top_n = 20)
#' #
#' # where the sign and magnitude of each gene's contribution are visualized.
#' }
plot_sign_importance <- function(fit = NULL, tab = NULL, top = 30,
                                 map_to_symbol = FALSE,
                                 from = "ENTREZID", to = "SYMBOL",
                                 show_legend = TRUE, palette = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plotting.", call. = FALSE)

  # Resolve the input table
  if (is.null(tab)) {
    if (is.null(fit))
      stop("Provide either `tab` (from sign_importance) or `fit` with signed_importance in fit@imp.", call. = FALSE)
    if (is.null(fit@imp) || !"signed_importance" %in% colnames(fit@imp))
      stop("`fit@imp` lacks `signed_importance`. Call sign_importance() first or pass `tab=`.", call. = FALSE)
    tab <- fit@imp
  }

  req_cols <- c("gene", "importance", "direction", "signed_importance")
  miss <- setdiff(req_cols, colnames(tab))
  if (length(miss))
    stop("`tab` is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)

  df <- tab[, req_cols, drop = FALSE]
  df <- df[order(abs(df$signed_importance), decreasing = TRUE), , drop = FALSE]
  df <- head(df, top)

  df$gene <- factor(df$gene, levels = rev(df$gene))
  df$direction <- as.character(df$direction)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene, y = signed_importance, fill = direction)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = sprintf("Top %d signed predictive genes", top),
      x = "Gene",
      y = "Signed importance"
    )

  if (!isTRUE(show_legend)) p <- p + ggplot2::theme(legend.position = "none")
  if (!is.null(palette))    p <- p + ggplot2::scale_fill_manual(values = palette)

  # Optional: map IDs -> SYMBOLs on x-axis like plot_importance()
  if (isTRUE(map_to_symbol)) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      warning("org.Hs.eg.db is not installed; showing original IDs.", call. = FALSE)
      return(p)
    }
    if (!"top_genes" %in% getNamespaceExports("RFGeneRank")) {
      warning("top_genes() not available; showing original IDs.", call. = FALSE)
      return(p)
    }
    tg <- RFGeneRank::top_genes(fit, n = top, map = TRUE, from = from, to = to)$table
    labs <- setNames(tg$mapped, tg$gene)
    p <- p + ggplot2::scale_x_discrete(labels = labs)
  }

  p
}

#' Expression-space embedding (PCA or UMAP) of top RF genes
#'
#' Plots samples in the original expression space using the top RF-ranked genes.
#' For PCA, axis labels include the percent variance explained.
#' For UMAP, the title shows the engine and key parameters used.
#'
#' @param fit  GeneRankFit
#' @param se   SummarizedExperiment with assay "expr"
#' @param n_top integer; number of top RF genes to use (default 100)
#' @param type  "umap" or "pca"
#' @param engine UMAP engine, "umap" or "uwot" (used when type="umap")
#' @param neighbors,min_dist,metric UMAP parameters
#' @param zscore Logical; z-score samples x genes matrix before embedding (default TRUE)
#' @param seed RNG seed for UMAP reproducibility (default \code{fit@params$seed})
#' @param point_size numeric
#' @param show_legend logical
#' @param palette optional named vector of colors (names must match levels)
#' @return A ggplot object.
#' @export
#' @importFrom SummarizedExperiment assay colData assayNames
#' @importFrom ggplot2 ggplot aes geom_point labs theme_classic scale_color_manual
#' @importFrom stats prcomp
#' @examples
#' # For reproducibility, specify a fixed seed (e.g., set.seed(1)) before running this example.
#'
#' # Toy expression: 15 genes × 8 samples
#' expr <- matrix(stats::rnorm(15 * 8), nrow = 15)
#' rownames(expr) <- paste0("gene", 1:15)
#' colnames(expr) <- paste0("sample", 1:8)
#'
#' # Binary labels for samples
#' label <- factor(rep(c("A", "B"), each = 4))
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(expr = expr),
#'   colData = data.frame(label = label)
#' )
#'
#' se
#'
#' \donttest{
#' # In practice, plot_embed_expr() uses the expression assay directly,
#' # for example:
#' #
#' #   plot_embed_expr(
#' #     se,
#' #     label_col = "label",
#' #     method    = "PCA"
#' #   )
#' #
#' # to visualize sample-level embeddings coloured by the outcome.
#' }
plot_embed_expr <- function(fit, se, n_top = 100,
                            type = c("umap","pca"),
                            engine = c("umap","uwot"),
                            neighbors = 15, min_dist = 0.1, metric = "euclidean",
                            zscore = TRUE, seed = fit@params$seed,
                            point_size = 2, show_legend = TRUE, palette = NULL) {
  type   <- match.arg(type)
  engine <- match.arg(engine)

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plotting.", call. = FALSE)
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("SummarizedExperiment is required.", call. = FALSE)

  X <- SummarizedExperiment::assay(se, "expr")
  md <- as.data.frame(SummarizedExperiment::colData(se))
  y  <- factor(md[[fit@params$label_col]])
  names(y) <- rownames(md)

  feats <- head(fit@imp$gene, n_top)
  feats <- intersect(feats, rownames(X))
  if (length(feats) < 2L)
    stop("Fewer than 2 overlapping top genes for embedding.", call. = FALSE)

  # samples x features
  M <- t(as.matrix(X[feats, , drop = FALSE]))
  if (isTRUE(zscore)) M <- scale(M)

  if (type == "pca") {
    # drop zero-variance columns for robustness
    if (!requireNamespace("matrixStats", quietly = TRUE))
      stop("matrixStats is required for PCA pre-check.", call. = FALSE)
    keep <- matrixStats::colSds(M) > 0
    if (sum(keep) < 2L) stop("Too few variable genes for PCA.", call. = FALSE)
    Mp <- M[, keep, drop = FALSE]

    pc   <- stats::prcomp(Mp, center = TRUE, scale. = FALSE)
    ve   <- (pc$sdev^2) / sum(pc$sdev^2)
    lab1 <- sprintf("PC1 (%.1f%%)", 100 * ve[1])
    lab2 <- sprintf("PC2 (%.1f%%)", 100 * ve[2])

    emb <- pc$x[, seq_len(2), drop = FALSE]
    df <- data.frame(PC1 = emb[, 1], PC2 = emb[, 2],
                     state = y[match(rownames(emb), names(y))],
                     sample = rownames(emb), row.names = NULL)

    p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = state)) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::theme_classic() +
      ggplot2::labs(title = sprintf("PCA of Top %d RF Genes", length(feats)),
                    x = lab1, y = lab2)

    if (!isTRUE(show_legend)) p <- p + ggplot2::theme(legend.position = "none")
    if (!is.null(palette))    p <- p + ggplot2::scale_color_manual(values = palette)
    return(p)
  }

  # ---- UMAP path ----
  if (engine == "umap") {
    if (!requireNamespace("umap", quietly = TRUE))
      stop("The 'umap' package is required for engine = 'umap'.", call. = FALSE)
    cfg <- umap::umap.defaults
    cfg$n_neighbors  <- neighbors
    cfg$min_dist     <- min_dist
    cfg$metric       <- metric
    cfg$random_state <- seed
    um  <- umap::umap(M, config = cfg)
    emb <- um$layout
    rownames(emb) <- rownames(M)
  } else {
    if (!requireNamespace("uwot", quietly = TRUE))
      stop("The 'uwot' package is required for engine = 'uwot'.", call. = FALSE)
    emb <- uwot::umap(M, n_neighbors = neighbors, min_dist = min_dist, metric = metric)
    rownames(emb) <- rownames(M)
  }

  df <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                   state = y[match(rownames(M), names(y))],
                   sample = rownames(M), row.names = NULL)

  p <- ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, color = state)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = sprintf("UMAP of Top %d RF Genes (%s, n=%d, min_dist=%.2f, metric=%s)",
                      length(feats), engine, neighbors, min_dist, metric),
      x = "UMAP1", y = "UMAP2"
    )

  if (!isTRUE(show_legend)) p <- p + ggplot2::theme(legend.position = "none")
  if (!is.null(palette))    p <- p + ggplot2::scale_color_manual(values = palette)
  p
}

#' Decision-space embedding (OOF probability space): PCA or UMAP
#'
#' Embeds samples using the out-of-fold (OOF) probability space (2D).
#' Good for inspecting decision geometry; for biology structure, prefer \code{plot_embed_expr()}.
#'
#' @param fit GeneRankFit (must have OOF probabilities).
#' @param type \code{"pca"} or \code{"umap"}.
#' @param engine UMAP engine if type="umap" (\code{"umap"} or \code{"uwot"}).
#' @param neighbors,min_dist,metric,seed UMAP params.
#' @param point_size numeric
#' @param show_legend logical
#' @param palette optional manual color palette
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point labs theme_classic scale_color_manual
#' @importFrom stats prcomp
#' @examples
#' # For reproducibility, specify a fixed seed (e.g., set.seed(1)) before running this example.
#'
#' # Toy 2D embedding for 10 samples
#' embed_df <- data.frame(
#'   sample = paste0("sample", 1:10),
#'   dim1   = stats::rnorm(10),
#'   dim2   = stats::rnorm(10),
#'   label  = factor(rep(c("A", "B"), each = 5))
#' )
#'
#' head(embed_df)
#'
#' \donttest{
#' # In practice, plot_embed() is used on a GeneRankFit object
#' # produced by the RFGeneRank workflow, for example:
#' #
#' #   fit <- gene_rank(se, genes = rownames(se), ...)
#' #   plot_embed(fit, palette = c("A" = "#1f77b4", "B" = "#ff7f0e"))
#' #
#' # where the embedding (e.g. PCA / UMAP) is stored inside 'fit'
#' # and coloured by the outcome or another covariate.
#' }
plot_embed <- function(fit,
                       type = c("pca","umap"),
                       engine = c("umap","uwot"),
                       neighbors = fit@params$umap$neighbors %||% 15,
                       min_dist = fit@params$umap$min_dist %||% 0.1,
                       metric   = fit@params$umap$metric   %||% "cosine",
                       seed     = fit@params$seed %||% 1,
                       point_size = 2, show_legend = TRUE, palette = NULL) {
  type   <- match.arg(type)
  engine <- match.arg(engine)

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plotting.", call. = FALSE)

  prob <- fit@oof$prob
  if (is.null(prob) || !is.matrix(prob))
    stop("fit@oof$prob is missing; run rank_genes() first.", call. = FALSE)

  # samples x classes matrix
  M <- prob
  y <- fit@oof$y

  if (type == "pca") {
    pc   <- stats::prcomp(M, center = TRUE, scale. = TRUE)
    ve   <- (pc$sdev^2) / sum(pc$sdev^2)
    lab1 <- sprintf("PC1 (%.1f%%)", 100 * ve[1])
    lab2 <- sprintf("PC2 (%.1f%%)", 100 * ve[2])

    df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                     state = y[match(rownames(M), names(y))])
    p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = state)) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::theme_classic() +
      ggplot2::labs(title = "Decision-space PCA (OOF probabilities)",
                    x = lab1, y = lab2)
  } else {
    if (engine == "umap") {
      if (!requireNamespace("umap", quietly = TRUE))
        stop("The 'umap' package is required for engine = 'umap'.", call. = FALSE)
      cfg <- umap::umap.defaults
      cfg$n_neighbors  <- neighbors
      cfg$min_dist     <- min_dist
      cfg$metric       <- metric
      cfg$random_state <- seed
      um  <- umap::umap(M, config = cfg)
      emb <- um$layout
      rownames(emb) <- rownames(M)
    } else {
      if (!requireNamespace("uwot", quietly = TRUE))
        stop("The 'uwot' package is required for engine = 'uwot'.", call. = FALSE)
      emb <- uwot::umap(M, n_neighbors = neighbors, min_dist = min_dist, metric = metric)
      rownames(emb) <- rownames(M)
    }
    df <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                     state = y[match(rownames(M), names(y))])
    p <- ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, color = state)) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title = sprintf("Decision-space UMAP (OOF) (%s, n=%d, min_dist=%.2f, metric=%s)",
                        engine, neighbors, min_dist, metric),
        x = "UMAP1", y = "UMAP2"
      )
  }

  if (!isTRUE(show_legend)) p <- p + ggplot2::theme(legend.position = "none")
  if (!is.null(palette))    p <- p + ggplot2::scale_color_manual(values = palette)
  p
}

# --- helpers used by plot_shap_dependence ------------------------------------

# coerce an age column to numeric safely
.make_age_numeric <- function(cd, age_col) {
  if (!age_col %in% names(cd)) stop("No '", age_col, "' found in colData.")
  x <- cd[[age_col]]
  if (is.numeric(x)) return(x)
  xnum <- try(as.numeric(as.character(x)), silent = TRUE)

  if (all(is.finite(xnum)))
    return(xnum)

  stop("Could not coerce '", age_col, "' to numeric.")
}

# equal-width breaks for age
.age_breaks <- function(x, width = 10, start = NULL) {
  if (is.null(start)) start <- floor(min(x, na.rm = TRUE) / width) * width
  seq(start, ceiling(max(x, na.rm = TRUE) / width) * width, by = width)
}

# parse "(a,b]" style labels from cut()
.parse_bin <- function(lbl) {
  s <- gsub("\\[|\\(|\\]|\\)", "", as.character(lbl))
  as.numeric(strsplit(s, ",")[[1]])
}

# summarize SHAP by age bins per group (for shading)
.summarize_age_bins <- function(df, width = 10, alpha = 0.05) {
  brks <- .age_breaks(df$x, width)
  tab <- df |>
    transform(age_bin = cut(x, breaks = brks, include.lowest = TRUE, right = FALSE)) |>
    aggregate(SHAP ~ age_bin + color,
              FUN = function(v) c(n = length(v), m = mean(v), sd = sd(v)))

  out <- do.call(data.frame, within(tab, {
    n      <- SHAP[, "n"];  mean_phi <- SHAP[, "m"];  sd_phi <- SHAP[, "sd"]
    se_phi <- sd_phi / sqrt(pmax(n, 1))
    t_stat <- ifelse(se_phi > 0 & n >= 2, mean_phi / se_phi, NA_real_)
    p_val  <- ifelse(is.finite(t_stat), 2 * stats::pt(-abs(t_stat), df = n - 1), NA_real_)
    fdr    <- stats::p.adjust(p_val, "BH")
    abs_t  <- abs(t_stat); abs_mean <- abs(mean_phi)
    SHAP   <- NULL
  }))
  attr(out, "alpha") <- alpha
  out
}

# add shaded Age x Sex bands (significant or top-k)
.add_band_highlight <- function(p, bins_tbl, alpha = 0.12,
                                mode = c("significant","top_k"), top_k = 1) {
  mode <- match.arg(mode)
  if (!nrow(bins_tbl)) return(p)

  pick <- if (mode == "significant") {
    a <- attr(bins_tbl, "alpha"); if (is.null(a)) a <- 0.05
    subset(bins_tbl, !is.na(fdr) & fdr <= a)
  } else {
    do.call(rbind, by(bins_tbl, bins_tbl$color, function(d) {
      d[order(-d$abs_t), ][seq_len(min(top_k, nrow(d))), , drop = FALSE]
    }))
  }
  if (!nrow(pick)) return(p)

  for (i in seq_len(nrow(pick))) {
    rng <- .parse_bin(pick$age_bin[i])
    fill_col <- if (as.character(pick$color[i]) %in% c("Female","female")) "#ff5ca8" else "#00bcd4"
    p <- p +
      ggplot2::annotate("rect", xmin = rng[1], xmax = rng[2], ymin = -Inf, ymax = Inf,
                        fill = fill_col, alpha = alpha, color = NA) +
      ggplot2::annotate("segment", x = mean(rng), xend = mean(rng), y = -Inf, yend = Inf,
                        linewidth = 0.25, alpha = 0.25)
  }
  p
}

# minimal palette used by the plot
.rfgr_sex_colors <- c(Female="#ff5ca8", Male="#00bcd4",
                      female="#ff5ca8", male="#00bcd4",
                      unknown="#7aa6ff", all="grey40")

# ---------- main: SHAP dependence scatter for one gene (cached-first, robust) ----------
#' SHAP dependence scatter for one gene (cached-first, robust, fast)
#'
#' @param fit GeneRankFit (must have fit@oof$prob and labels in colData(se)$state).
#' @param se  SummarizedExperiment with assay "expr".
#' @param gene Character gene ID present in assay(se,"expr").
#' @param x    Covariate name in colData(se) for x-axis (default "age").
#' @param color_by Optional grouping column name in colData(se); if NULL, no colour grouping.
#' @param palette Optional character vector of colour names or hex codes.
#'   If \code{color_by} is:
#'   \itemize{
#'     \item numeric: \code{palette} (length >= 2) is used as a continuous gradient via \code{scale_color_gradientn()}.
#'     \item categorical: \code{palette} is recycled/trimmed to the number of levels and used via \code{scale_color_manual()}.
#'   }
#'   If \code{palette} is NULL, ggplot2 defaults are used, except for \code{color_by = "sex"} /
#'   \code{"sex_clean"} where a fixed pink/cyan palette is used.
#' @param shap_mat Optional SHAP matrix (rows = samples, cols = genes).
#' @param nsim Fastshap nsim if auto-computing (default 64).
#' @param pos_label Positive class label (defaults to level 2 of outcome).
#' @param age_band_width Width of shaded age bands (years) for sex plots.
#' @param band_alpha Alpha of bands.
#' @param band_mode "significant" (FDR<=0.05) or "top_k".
#' @param top_k If band_mode="top_k", how many bins per group to shade.
#'
#' @return ggplot object.
#' @export
#' @examples
#' expr <- matrix(
#'   rnorm(10 * 20),
#'   nrow = 10,
#'   dimnames = list(
#'     paste0("gene", 1:10),
#'     paste0("sample", 1:20)
#'   )
#' )
#'
#' y <- factor(rep(c("Control", "Case"), each = 10))
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(expr = expr),
#'   colData = data.frame(
#'     state = y,
#'     age   = seq_len(20) + 40
#'   )
#' )
#'
#' fit <- rank_genes(se, label_col = "state", trees = 50)
#'
#' # Minimal SHAP-like matrix: one column named by the gene
#' shap_mat <- matrix(
#'   rnorm(20),
#'   nrow = 20,
#'   dimnames = list(NULL, "gene1")
#' )
#'
#' plot_shap_dependence(
#'   fit      = fit,
#'   se       = se,
#'   gene     = "gene1",
#'   x        = "age",
#'   shap_mat = shap_mat
#' )
plot_shap_dependence <- function(fit, se,
                                 gene,
                                 x = "age",
                                 color_by = NULL,
                                 palette = NULL,
                                 shap_mat = NULL,
                                 nsim = 64,
                                 pos_label = NULL,
                                 age_band_width = 10,
                                 band_alpha = 0.12,
                                 band_mode = c("significant","top_k"),
                                 top_k = 1) {

  band_mode <- match.arg(band_mode)
  stopifnot("expr" %in% SummarizedExperiment::assayNames(se))

  ## ---- align to OOF order ---------------------------------------------------
  sids <- rownames(fit@oof$prob)
  cd   <- as.data.frame(SummarizedExperiment::colData(se))[sids, , drop = FALSE]
  yfac <- droplevels(as.factor(cd$state))
  if (is.null(pos_label)) pos_label <- levels(yfac)[2]

  if (!x %in% names(cd))
    stop(sprintf("covariate '%s' not found in colData(se)", x))
  if (!is.null(color_by) && !color_by %in% names(cd))
    stop(sprintf("grouping covariate '%s' not found in colData(se)", color_by))

  ## ---- SHAP (prefer precomputed matrix) -------------------------------------
  if (!is.null(shap_mat)) {
    if (is.null(rownames(shap_mat)) && nrow(shap_mat) == length(sids))
      rownames(shap_mat) <- sids
    common <- intersect(sids, rownames(shap_mat))
    if (!length(common))
      stop("No overlap between shap_mat rownames and OOF IDs.")
    shap_mat <- shap_mat[common, , drop = FALSE]
    cd       <- cd[common, , drop = FALSE]
  } else {
    if (!requireNamespace("fastshap", quietly = TRUE))
      stop("Provide `shap_mat=` or install 'fastshap' to auto-compute SHAP.", call. = FALSE)

    X <- t(SummarizedExperiment::assay(se, "expr"))[sids, , drop = FALSE]
    if (!gene %in% colnames(X))
      stop(sprintf("gene '%s' not found in assay(se,'expr').", gene))

    pred_fun <- function(object, newdata) {
      nd <- as.data.frame(newdata)

      pr <- tryCatch(stats::predict(object, newdata = nd, type = "prob"),
                     error = function(e) NULL)
      if (!is.null(pr)) return(as.numeric(pr[, pos_label, drop = TRUE]))

      pr <- tryCatch(stats::predict(object, newdata = nd, type = "response"),
                     error = function(e) NULL)
      if (!is.null(pr)) {
        if (is.matrix(pr) || is.data.frame(pr)) return(as.numeric(pr[, pos_label, drop = TRUE]))
        return(as.numeric(pr))
      }

      pr <- tryCatch(predict(object, data = nd, type = "response")$predictions,
                     error = function(e) NULL)
      if (!is.null(pr)) {
        if (is.matrix(pr) || is.data.frame(pr)) return(as.numeric(pr[, pos_label, drop = TRUE]))
        return(as.numeric(pr))
      }

      stop("pred_wrapper: could not obtain positive-class probabilities.")
    }

    X_full <- as.data.frame(X)
    phi <- fastshap::explain(
      object        = fit,
      X             = X_full,
      newdata       = X_full,
      nsim          = nsim,
      feature_names = gene,
      pred_wrapper  = pred_fun
    )
    shap_mat <- as.matrix(phi)
    if (!gene %in% colnames(shap_mat)) colnames(shap_mat) <- gene
    rownames(shap_mat) <- rownames(X)
  }

  if (!gene %in% colnames(shap_mat))
    stop(sprintf("Gene '%s' not present in SHAP matrix.", gene))

  ## ---- build plotting data.frame -------------------------------------------
  df <- data.frame(
    x    = cd[[x]],
    SHAP = shap_mat[, gene],
    row.names = rownames(cd)
  )

  if (is.null(color_by)) {
    df$group <- factor("all")
  } else {
    df$group <- cd[[color_by]]
  }

  ## ---- age coercion if requested -------------------------------------------
  x_is_age <- x %in% c("age", "age_num") || grepl("^age$", x)
  if (x_is_age)
    df$x <- .make_age_numeric(cd, x)

  ## ---- clean NA / non-finite -----------------------------------------------
  df <- df[is.finite(df$x) & is.finite(df$SHAP), , drop = FALSE]
  has_group <- !is.null(color_by)
  if (has_group)
    df <- df[!is.na(df$group), , drop = FALSE]

  ## ---- grouping metadata ----------------------------------------------------
  sex_like    <- has_group && tolower(color_by) %in% c("sex","sex_clean")
  group_is_num <- has_group && is.numeric(df$group)

  ## ---- base geometry --------------------------------------------------------
  if (is.numeric(df$x)) {

    show_legend_points <- !(sex_like && x_is_age)

    if (has_group) {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = SHAP, color = group)) +
        ggplot2::geom_point(alpha = 0.9, size = 2.0, show.legend = show_legend_points)
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = SHAP)) +
        ggplot2::geom_point(alpha = 0.9, size = 2.0, color = "grey40")
    }

    p <- p +
      ggplot2::labs(
        x = x,
        y = sprintf("SHAP contribution of %s to P(%s)", gene, pos_label)
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin      = ggplot2::margin(5, 5, 5, 5)
      )

    ## colour scales based on palette -----------------------------------------
    if (has_group) {
      if (group_is_num) {
        if (!is.null(palette)) {
          if (length(palette) < 2)
            stop("For numeric 'color_by', 'palette' must have at least 2 colours.")
          p <- p + ggplot2::scale_color_gradientn(colours = palette)
        } else {
          p <- p + ggplot2::scale_color_gradient(low = "#FFA726", high = "#B71C1C")
        }
      } else {
        df$group <- droplevels(factor(df$group))
        levs <- levels(df$group)

        if (sex_like && is.null(palette)) {
          pal_use <- .rfgr_sex_colors[levs]
          pal_use[is.na(pal_use)] <- "grey70"
        } else if (!is.null(palette)) {
          pal_use <- rep(palette, length.out = length(levs))
          names(pal_use) <- levs
        } else {
          pal_use <- NULL
        }

        if (!is.null(pal_use)) {
          guide_mode <- if (sex_like && x_is_age) "none" else "legend"
          p <- p + ggplot2::scale_color_manual(values = pal_use, guide = guide_mode)
        }
      }
    }

    ## smoothing logic --------------------------------------------------------
    nux <- length(unique(df$x))
    if (nux >= 8 && requireNamespace("mgcv", quietly = TRUE)) {
      p <- p + ggplot2::geom_smooth(
        se      = FALSE,
        method  = mgcv::gam,
        formula = y ~ s(x, k = 6)
      )
    } else if (nux >= 8) {
      p <- p + ggplot2::geom_smooth(
        se      = FALSE,
        method  = "loess",
        formula = y ~ x,
        span    = 0.75
      )
    } else {
      df_sum <- aggregate(
        SHAP ~ group + x, data = df,
        FUN = function(v) c(mean = mean(v), se = stats::sd(v) / sqrt(length(v)))
      )
      df_sum$mean <- df_sum$SHAP[, "mean"]
      df_sum$se   <- df_sum$SHAP[, "se"]
      df_sum$SHAP <- NULL

      p <- p +
        ggplot2::geom_errorbar(
          data        = df_sum,
          ggplot2::aes(ymin = mean - se, ymax = mean + se),
          width       = 0,
          linewidth   = 0.5,
          inherit.aes = FALSE
        ) +
        ggplot2::geom_line(
          data        = df_sum,
          ggplot2::aes(y = mean, group = group, color = group),
          linewidth   = 1,
          inherit.aes = FALSE
        )
    }

  } else {
    ## categorical x: violin/boxplot form ------------------------------------
    df$x <- factor(df$x)
    if (has_group) {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = SHAP, fill = group)) +
        ggplot2::geom_violin(trim = TRUE) +
        ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.5)
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = SHAP)) +
        ggplot2::geom_violin(trim = TRUE, fill = "grey80") +
        ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.5)
    }

    p <- p +
      ggplot2::labs(
        x = x,
        y = sprintf("SHAP contribution of %s to P(%s)", gene, pos_label)
      ) +
      ggplot2::theme_classic(base_size = 13)

    if (has_group && !group_is_num) {
      df$group <- droplevels(factor(df$group))
      levs <- levels(df$group)

      if (sex_like && is.null(palette)) {
        pal_use <- .rfgr_sex_colors[levs]
        pal_use[is.na(pal_use)] <- "grey70"
      } else if (!is.null(palette)) {
        pal_use <- rep(palette, length.out = length(levs))
        names(pal_use) <- levs
      } else {
        pal_use <- NULL
      }

      if (!is.null(pal_use)) {
        p <- p + ggplot2::scale_fill_manual(values = pal_use)
      }
    }
  }

  ## ---- Age x Sex band highlighting (only for sex on age) ---------------------
  wants_bands <- is.numeric(df$x) &&
    x_is_age &&
    sex_like &&
    nrow(df) > 0 &&
    length(unique(df$group)) >= 2

  if (wants_bands) {
    df_bands <- data.frame(
      x     = df$x,
      SHAP  = df$SHAP,
      color = df$group
    )
    df_bands <- df_bands[is.finite(df_bands$x) & !is.na(df_bands$color), , drop = FALSE]

    if (nrow(df_bands) > 0) {
      bins_tbl <- .summarize_age_bins(
        df_bands,
        width = age_band_width,
        alpha = 0.05
      )
      p <- .add_band_highlight(
        p, bins_tbl,
        alpha = band_alpha,
        mode  = band_mode,
        top_k = top_k
      )
    }

    if (requireNamespace("patchwork", quietly = TRUE)) {
      key_df <- data.frame(
        group = factor(c("Female","Male"), levels = c("Female","Male")),
        y     = c(1.5, 0.5)
      )
      p_key <- ggplot2::ggplot(key_df, ggplot2::aes(1, y, fill = group)) +
        ggplot2::geom_tile(height = 1.0, width = 0.35) +
        ggplot2::scale_fill_manual(values = .rfgr_sex_colors, guide = "none") +
        ggplot2::geom_text(
          ggplot2::aes(label = group),
          angle = 90,
          vjust = 0.5,
          hjust = 0.5,
          size  = 4.2
        ) +
        ggplot2::coord_cartesian(xlim = c(0.8, 2.2), ylim = c(0, 2), clip = "off") +
        ggplot2::labs(y = "Sex", x = NULL) +
        ggplot2::theme_void(base_size = 14) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_text(
            angle  = 90,
            vjust  = 0.5,
            margin = ggplot2::margin(r = 4)
          ),
          plot.margin = ggplot2::margin(5, 10, 5, 0)
        )

      return(p + p_key + patchwork::plot_layout(widths = c(1, 0.12)))
    }
  }

  p
}

#' Compute a SHAP matrix using a single ranger model
#'
#' Trains a single probability random forest on the selected genes and computes
#' per-sample SHAP values with fastshap. Rows = samples (colData rows), cols = genes.
#'
#' @param se SummarizedExperiment with assay "expr" (genes x samples).
#' @param label_col Outcome column in colData(se) (factor).
#' @param genes Character vector of gene IDs (must match rownames of assay(se,"expr")).
#' @param class_weights Optional named numeric vector of class weights; if NULL,
#'   inverse-frequency weights are used.
#' @param num.trees Number of trees (default 500).
#' @param seed RNG seed.
#' @param nsim fastshap Monte Carlo samples (default 64).
#' @param pos_label Positive class label (defaults to level 2 of y).
#' @return Numeric matrix of SHAP values (n_samples x length(genes)).
#' @export
#' @importFrom SummarizedExperiment assay colData assayNames
#' @examples
#' # For reproducibility, specify a fixed seed (e.g., set.seed(1)) before running this example.
#'
#' # Toy feature matrix and binary outcome
#' X <- data.frame(
#'   x1 = stats::rnorm(20),
#'   x2 = stats::rnorm(20)
#' )
#' y <- factor(rep(c("A", "B"), each = 10))
#'
#' # Inspect inputs
#' head(X)
#'
#' \donttest{
#' # In practice, shap_train_ranger() computes a SHAP matrix using a single
#' # ranger model. A typical call is:
#' #
#' #   shap_obj <- shap_train_ranger(
#' #     X    = X,
#' #     y    = y,
#' #     nsim = 128
#' #   )
#' #
#' # where 'shap_obj' contains the SHAP values for each feature.
#' }
shap_train_ranger <- function(se, label_col = "state", genes,
                              class_weights = NULL,
                              num.trees = 500, seed = 1L, nsim = 64,
                              pos_label = NULL) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("Package 'ranger' is required.", call. = FALSE)
  if (!requireNamespace("fastshap", quietly = TRUE))
    stop("Package 'fastshap' is required.", call. = FALSE)

  stopifnot("expr" %in% SummarizedExperiment::assayNames(se))
  Xall <- t(as.matrix(SummarizedExperiment::assay(se, "expr")))  # samples x genes
  y    <- droplevels(as.factor(SummarizedExperiment::colData(se)[[label_col]]))
  if (is.null(pos_label)) pos_label <- levels(y)[2]

  genes <- intersect(genes, colnames(Xall))
  if (length(genes) < 1)
    stop("None of the requested genes are present in assay(se,'expr').")

  X <- Xall[, genes, drop = FALSE]
  dat <- data.frame(y = y, X, check.names = FALSE)

  # class weights (inverse frequency by default)
  if (is.null(class_weights)) {
    tab <- table(y)
    class_weights <- as.numeric(median(tab) / tab)
    names(class_weights) <- names(tab)
  }
  wts <- class_weights[as.character(y)]

  fit_rg <- ranger::ranger(
    formula     = y ~ .,
    data        = dat,
    probability = TRUE,
    num.trees   = num.trees,
    case.weights = wts
  )

  # fastshap predictor: return probability for positive class
  pred_fun <- function(object, newdata) {
    pr <- predict(object, data = newdata, type = "response")$predictions
    if (is.matrix(pr)) pr <- pr[, pos_label, drop = TRUE]
    as.numeric(pr)
  }

  phi <- fastshap::explain(
    object        = fit_rg,
    X             = as.data.frame(X),
    newdata       = as.data.frame(X),
    nsim          = nsim,
    feature_names = colnames(X),
    pred_wrapper  = pred_fun
  )

  SHAP <- as.matrix(phi)
  rownames(SHAP) <- rownames(X)
  SHAP
}

#' ROC curve from OOF probabilities (single model)
#' @param fit GeneRankFit with \code{fit@oof$prob} and \code{fit@oof$y}
#' @return A \code{ggplot2} object containing the ROC curve derived from
#' @export
#' @importFrom pROC roc auc coords
#' @importFrom ggplot2 ggplot aes geom_step geom_abline labs theme_bw annotate
#' @examples
#' expr <- matrix(
#'   rnorm(10 * 20),
#'   nrow = 10,
#'   dimnames = list(
#'     paste0("gene", 1:10),
#'     paste0("sample", 1:20)
#'   )
#' )
#'
#' y <- factor(rep(c("Control", "Case"), each = 10))
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(expr = expr),
#'   colData = data.frame(state = y)
#' )
#'
#' fit <- rank_genes(se, label_col = "state", trees = 100)
#' plot_roc(fit)
plot_roc <- function(fit) {
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("pROC is required for ROC plotting.", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plotting.", call. = FALSE)

  prob <- fit@oof$prob
  y    <- droplevels(as.factor(fit@oof$y))
  stopifnot(is.matrix(prob), nrow(prob) == length(y))

  # Positive class = last level
  pos <- tail(levels(y), 1)
  p   <- prob[, pos]

  roc <- pROC::roc(response = y, predictor = p,
                   levels = c(setdiff(levels(y), pos), pos),
                   direction = "<")

  crd <- pROC::coords(roc, "all", ret = c("specificity","sensitivity"), transpose = FALSE)
  df  <- data.frame(FPR = 1 - crd$specificity, TPR = crd$sensitivity)
  df  <- df[order(df$FPR, df$TPR), , drop = FALSE]
  df  <- df[!duplicated(df), , drop = FALSE]

  auc_val <- as.numeric(pROC::auc(roc))
  auc_txt <- sprintf("AUC = %.3f", auc_val)

  ggplot2::ggplot(df, ggplot2::aes(FPR, TPR)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", linewidth = 0.7, color = "gray55") +
    ggplot2::geom_step(direction = "vh",
                       linewidth = 1.6, lineend = "round", color = "red") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::scale_x_continuous(name = "1 - Specificity (FPR)", breaks = seq(0, 1, 0.25)) +
    ggplot2::scale_y_continuous(name = "Sensitivity (TPR)", breaks = seq(0, 1, 0.25)) +
    ggplot2::annotate("text", x = 0.65, y = 0.25, label = auc_txt, size = 5) +
    ggplot2::labs(title = "OOF ROC") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "gray90")
    )
}

# ==============================================================================
# 2) Auxiliary helper functions
# ==============================================================================

#' Multi-model ROC (guardrailed OOF)
#'
#' Overlays ROC curves for one or more fits that expose \code{$oof} with
#' \code{p_use} and \code{y}. Works with \code{validate_genes()} outputs.
#'
#' @param fits Named list of model results (e.g., \code{list(ranger=val$ranger, ...)}).
#' @param title Character plot title.
#' @return A ggplot object.
#' @export
#' @importFrom pROC roc auc coords
#' @importFrom ggplot2 ggplot aes geom_abline geom_step labs theme_minimal scale_x_continuous scale_y_continuous
#' @examples
#' # Toy binary outcome
#' y <- factor(rep(c("Control", "Case"), each = 10))
#'
#' # Positive-class probabilities (Case = positive class)
#' p_use1 <- runif(20)
#' p_use2 <- runif(20)
#'
#' # Create minimal fit objects expected by plot_roc_multi()
#' fits <- list(
#'   Model_1 = list(
#'     oof = list(
#'       y = y,
#'       p_use = p_use1
#'     )
#'   ),
#'   Model_2 = list(
#'     oof = list(
#'       y = y,
#'       p_use = p_use2
#'     )
#'   )
#' )
#'
#' plot_roc_multi(fits)
plot_roc_multi <- function(fits, title = "ROC (OOF, guardrailed)") {
  if (!requireNamespace("pROC", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE))
    stop("Packages 'pROC' and 'ggplot2' are required.", call. = FALSE)

  mk_roc_df <- function(x, name) {
    oof <- x$oof
    if (is.null(oof$p_use)) stop("oof$p_use missing; run validate_genes() first.")
    lev <- levels(oof$y)
    roc <- pROC::roc(oof$y, oof$p_use, levels = rev(lev), quiet = TRUE)
    crd <- pROC::coords(roc, "all", ret = c("specificity","sensitivity"), transpose = FALSE)
    data.frame(
      fpr = 1 - crd$specificity,
      tpr = crd$sensitivity,
      method = name,
      auc = as.numeric(pROC::auc(roc)),
      stringsAsFactors = FALSE
    )
  }
  dfl <- Map(mk_roc_df, fits, names(fits))
  df  <- do.call(rbind, dfl)

  ggplot2::ggplot(df, ggplot2::aes(fpr, tpr, color = method)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::geom_step(linewidth = 1, direction = "vh") +
    ggplot2::labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = title, color = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
    ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25))
}

#' Confusion-matrix heatmap
#'
#' @param cm A 2x2 table as returned by \code{val$<method>$conf_mat}.
#' @param mode Either \code{"counts"} or \code{"rowpct"}.
#' @return A ggplot object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile geom_text coord_equal labs theme_minimal scale_fill_gradient
#' @importFrom scales percent
#' @examples
#' # For reproducibility, specify a fixed seed (e.g., set.seed(1)) before running this example.
#'
#' true <- factor(rep(c("A", "B"), each = 5))
#' pred <- factor(sample(c("A", "B"), 10, replace = TRUE))
#'
#' tab <- table(True = true, Predicted = pred)
#'
#' plot_confusion_heatmap(tab)
plot_confusion_heatmap <- function(cm, mode = c("counts","rowpct")) {
  mode <- match.arg(mode)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required.", call. = FALSE)

  df <- as.data.frame(cm)
  names(df) <- c("True","Predicted","Freq")
  if (mode == "rowpct") {
    df$pct <- with(df, Freq / ave(Freq, True, FUN = sum))
    lab <- ggplot2::aes(label = scales::percent(pct, accuracy = 0.1))
    fill <- ggplot2::aes(fill = pct)
    scale <- ggplot2::scale_fill_gradient(labels = scales::percent, low = "steelblue", high = "darkred")
    title <- "Confusion Matrix (row %)"
  } else {
    lab <- ggplot2::aes(label = Freq)
    fill <- ggplot2::aes(fill = Freq)
    scale <- ggplot2::scale_fill_gradient(low = "steelblue", high = "darkred")
    title <- "Confusion Matrix (counts)"
  }

  ggplot2::ggplot(df, ggplot2::aes(True, Predicted)) +
    ggplot2::geom_tile(fill) +
    ggplot2::geom_text(lab, color = "white", size = 5) +
    scale +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = "True", y = "Predicted") +
    ggplot2::theme_minimal(base_size = 12)
}

#' One-call plot suite (optional file export)
#'
#' Produces: importance, signed-importance, expression PCA, decision PCA,
#' ROC overlay (if \code{val} provided), confusion heatmaps, and optional SHAP.
#'
#' @param fit GeneRankFit from \code{rank_genes()}.
#' @param se  SummarizedExperiment with assay "expr".
#' @param val Optional result from \code{validate_genes()}.
#' @param outdir Optional directory to save PNGs (if not NULL).
#' @param top Integer; number of genes in importance plots.
#' @param shap_gene Optional gene ID for SHAP dependence.
#' @return (Invisibly) list of ggplot objects.
#' @export
#' @importFrom ggplot2 ggsave
#' @examples
#' # Toy expression matrix: genes x samples
#' expr <- matrix(
#'   rnorm(10 * 12),
#'   nrow = 10
#' )
#'
#' # Use known human Entrez IDs as gene names so annotation works
#' rownames(expr) <- c(
#'   "1",    # A1BG
#'   "2",    # A2M
#'   "9",    # NAT1
#'   "1956", # EGFR
#'   "2064", # ERBB2
#'   "5290", # PIK3CA
#'   "5728", # PTEN
#'   "7422", # VEGFA
#'   "1950", # EDN1
#'   "7157"  # TP53
#' )
#'
#' colnames(expr) <- paste0("sample", 1:12)
#'
#' # Binary phenotype stored in 'label' column
#' label <- factor(rep(c("Control", "Case"), each = 6))
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(expr = expr),
#'   colData = data.frame(
#'     label     = label,
#'     row.names = colnames(expr)
#'   )
#' )
#'
#' # Fit a small GeneRank model
#' fit <- rank_genes(
#'   se        = se,
#'   label_col = "label",
#'   n_top     = 10,
#'   trees     = 100
#' )
#'
#' # For this example, ensure fit@imp has signed_importance information
#' # required by plot_sign_importance().
#' if (is.null(fit@imp)) {
#'   ng <- nrow(expr)
#'   fit@imp <- data.frame(
#'     gene             = rownames(expr),
#'     importance       = seq_len(ng),
#'     direction        = rep(1L, ng),
#'     signed_importance = seq_len(ng),
#'     stringsAsFactors = FALSE
#'   )
#' } else {
#'   if (!"direction" %in% colnames(fit@imp)) {
#'     fit@imp$direction <- 1L
#'   }
#'   if (!"signed_importance" %in% colnames(fit@imp)) {
#'     fit@imp$signed_importance <- fit@imp$importance * fit@imp$direction
#'   }
#' }
#'
#' # Generate a suite of diagnostic plots (returned as a list of ggplot objects)
#' plots <- rfgr_plot_suite(fit, se)
#'
#' # Inspect available plots
#' names(plots)
rfgr_plot_suite <- function(fit, se, val = NULL, outdir = NULL, top = 30, shap_gene = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required.", call. = FALSE)

  plots <- list()

  # 1) importance / signed-importance
  p_imp <- plot_importance(fit, top = top, map_to_symbol = TRUE)
  p_sig <- plot_sign_importance(fit, top = top, map_to_symbol = TRUE)
  plots$importance <- p_imp
  plots$signed_importance <- p_sig

  # 2) expression-space embedding (top RF genes): PCA by default
  p_pca_expr <- plot_embed_expr(fit, se, n_top = 100, type = "pca")
  plots$expr_pca <- p_pca_expr

  # 3) decision-space embedding from OOF probabilities
  p_dec <- plot_embed(fit, type = "pca")
  plots$decision_pca <- p_dec

  # 4) multi-ROC from validate_genes() (if provided) + confusion heatmaps
  if (!is.null(val)) {
    fits <- setNames(lapply(val$summary$method, function(m) val[[m]]), val$summary$method)
    plots$roc_multi <- plot_roc_multi(fits, title = "ROC (OOF, guardrailed)")
    for (m in names(fits)) {
      plots[[paste0("cm_", m, "_counts")]]  <- plot_confusion_heatmap(fits[[m]]$conf_mat, "counts")
      plots[[paste0("cm_", m, "_rowpct")]]  <- plot_confusion_heatmap(fits[[m]]$conf_mat, "rowpct")
    }
  }

  # 5) optional SHAP dependence
  if (!is.null(shap_gene)) plots$shap <- plot_shap_dependence(fit, se, shap_gene)

  # optional save
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    save_plot <- function(p, nm, w=8, h=5) {
      fn <- file.path(outdir, paste0(nm, ".png"))
      ggplot2::ggsave(fn, p, dpi = 300, width = w, height = h)
    }
    save_plot(p_imp, "importance_top")
    save_plot(p_sig, "importance_signed_top")
    save_plot(p_pca_expr, "expr_pca")
    save_plot(p_dec, "decision_pca")
    if (!is.null(val)) {
      save_plot(plots$roc_multi, "roc_multi", w=7, h=5)
      for (nm in names(plots)) if (startsWith(nm, "cm_")) save_plot(plots[[nm]], nm, w=5, h=4)
    }
    if ("shap" %in% names(plots)) save_plot(plots$shap, "shap_dependence", w=7, h=4.5)
  }

  invisible(plots)
}
