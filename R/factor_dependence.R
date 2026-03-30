#' Covariate dependence of gene contributions (SHAP/proxy; "expr" assay)
#'
#' @param fit            Trained object containing out-of-fold predictions (n x k) and labels.
#'                       For S4 wrappers (e.g., GeneRankFit), the finalized learner should be in
#'                       A final fitted model and the training feature names may also be present.
#'                       this function will train a temporary full-data model on-the-fly (not saved).
#' @param se             SummarizedExperiment with assay "expr" (genes x samples).
#' @param covariates     Character vector of covariate names in colData(se), e.g. c("sex","age").
#' @param method         "shap" (default; uses fastshap if available) or "proxy".
#' @param ngenes         Number of genes to test (default 500). Use "ALL" or NULL for all genes.
#' @param gene_selection "importance" (default) or "variance".
#' @param nsim           SHAP Monte-Carlo permutations (default 128).
#' @param bg_per_class   Max background samples per class for SHAP (default 64).
#' @param cache_dir      Cache directory for SHAP matrices; default R user cache dir.
#' @param fdr_method     Multiple-testing correction (default "BH").
#' @param seed           RNG seed (default 1).
#' @param pred_fun       OPTIONAL user predictor: function(newdata) -> numeric p(positive).
#'                       If provided, it takes precedence over the built-in predictor.
#' @param pos_label      OPTIONAL (legacy) positive class label; kept for back-compat.
#' @param positive       OPTIONAL override for positive class. Either a class label (character)
#'                       or an index 1/2. If provided, it overrides `pos_label`.
#'
#' @return data.frame with columns: gene, covariate, test, stat, pval, fdr, effect, dependent.
#'         Attributes: "method","assay","genes_used","cache_file".
#' @importFrom stats aggregate sd var predict lm anova cor.test p.adjust resid
#' @examples
#' #' # For reproducibility, set a fixed seed such as set.seed(1) before running this example.
#'
#' # Tiny expression matrix
#' expr <- matrix(stats::rnorm(30), nrow = 6)
#' rownames(expr) <- paste0("gene", 1:6)
#' colnames(expr) <- paste0("sample", 1:5)
#'
#' # Covariates
#' cd <- data.frame(
#'   label = factor(c("A", "A", "B", "B", "B")),
#'   sex   = factor(c("M", "F", "F", "M", "M")),
#'   age   = c(30, 40, 35, 50, 60)
#' )
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(expr = expr),
#'   colData = cd
#' )
#'
#' # Minimal mock GeneRankFit object with required slots
#' prob_mat <- matrix(
#'   stats::runif(5),
#'   ncol = 1,
#'   dimnames = list(colnames(se), "B")
#' )
#'
#' fit <- methods::new(
#'   "GeneRankFit",
#'   oof = list(
#'     prob = prob_mat,
#'     y    = cd$label
#'   )
#' )
#'
#' @export
factor_dependence <- function(
  fit, se, covariates,
  method = c("shap","proxy"),
  ngenes = 500L,
  gene_selection = c("importance","variance"),
  nsim = 128L,
  bg_per_class = 64L,
  cache_dir = NULL,
  fdr_method = "BH",
  seed = 1L,
  pred_fun = NULL,
  pos_label = NULL,   # kept for compatibility
  positive = NULL     # new optional override
) {
  oof0 <- .rfgr_oof(fit)
  imp0 <- .rfgr_imp(fit)
  pars <- .rfgr_params(fit)
  mdl0 <- .rfgr_final_model(fit)
  fts0 <- .rfgr_features(fit)
  method <- match.arg(method)
  gene_selection <- match.arg(gene_selection)
  stopifnot(length(covariates) >= 1)

  ## ---------- Assay & labels ----------
  if (!"expr" %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay 'expr' not found in `se`. Ensure prepare_data() stored the matrix as 'expr'.")
  }
  expr_all <- SummarizedExperiment::assay(se, "expr")  # genes x samples
  cd_all   <- as.data.frame(SummarizedExperiment::colData(se))

  sids <- rownames(oof0$prob)
  if (is.null(sids)) stop("OOF probabilities must have rownames = sample IDs.")
  if (!all(sids %in% colnames(expr_all))) {
    miss <- setdiff(sids, colnames(expr_all))
    stop("Sample ID mismatch: some OOF sample IDs not in assay 'expr'. Examples: ",
         paste(utils::head(miss, 5), collapse = ", ")) # show only the first few missing items to keep the diagnostic message readable
  }
  if (!all(sids %in% rownames(cd_all))) {
    miss <- setdiff(sids, rownames(cd_all))
    stop("Sample ID mismatch: some OOF sample IDs not in colData(se). Examples: ",
         paste(utils::head(miss, 5), collapse = ", ")) # show only the first few missing items to keep the diagnostic message readable
  }

  expr <- expr_all[, sids, drop = FALSE]  # # genes x samples (aligned)
  cd   <- cd_all[sids, , drop = FALSE]
  if (!"state" %in% names(cd)) stop("colData(se)$state is required.")

  y <- droplevels(as.factor(cd$state))
  if (nlevels(y) != 2L) stop("Binary classification expected in colData(se)$state.")

  ## ---------- Positive class resolution ----------
  if (!is.null(positive)) {
    if (is.character(positive) && length(positive) == 1) {
      pos <- positive
    } else if (is.numeric(positive) && length(positive) == 1 && positive %in% c(1,2)) {
      pos <- levels(y)[positive]
    } else {
      stop("`positive` must be NULL, a single class label, or numeric 1/2.")
    }
  } else if (!is.null(pos_label)) {
    pos <- pos_label
  } else {
    pos <- levels(y)[2]  # default convention
  }
  if (!(pos %in% levels(y))) {
    stop("Chosen positive class ('", pos, "') is not in levels(y): ",
         paste(levels(y), collapse = ", "))
  }

  # don't allow the response label among covariates
  covariates <- setdiff(covariates, c("state","disease","label"))

  ## ---------- Gene set ----------
  Xg_full <- expr
  if (is.character(ngenes) && toupper(ngenes) == "ALL" || is.null(ngenes)) {
    genes <- rownames(Xg_full)
  } else if (gene_selection == "importance") {
    if (is.null(imp0) || !all(c("gene","importance") %in% names(imp0))) {
    stop("`gene_selection='importance'` requires an importance table with columns 'gene' and 'importance'.")
    }
    imp_df <- imp0[order(-imp0$importance), c("gene","importance"), drop = FALSE]
    genes  <- head(imp_df$gene, ngenes) # select the top-ranked genes for downstream dependence analysis
    genes  <- intersect(genes, rownames(Xg_full))
    if (!length(genes)) stop("No overlap between top-importance genes and 'expr' rownames.")
  } else { # variance
    v <- matrixStats::rowVars(as.matrix(Xg_full))
    names(v) <- rownames(Xg_full)
    ord   <- order(v, decreasing = TRUE)
    genes <- rownames(Xg_full)[ord][seq_len(min(ngenes, length(ord)))]
    genes <- genes[is.finite(v[genes]) & v[genes] > 0]
    if (!length(genes)) stop("No variable genes found in assay 'expr'.")
  }

  # Evaluation design (samples x genes)
  X_eval <- t(Xg_full[genes, , drop = FALSE])

  contr <- NULL
  used  <- "proxy"
  cache_path <- NA_character_

  ## ---------- SHAP path ----------
  if (method == "shap" && requireNamespace("fastshap", quietly = TRUE)) {
    used <- "shap"

    bg_idx <- unlist(lapply(split(seq_len(nrow(cd)), cd$state), function(ix) {
      sample(ix, size = min(bg_per_class, length(ix)))
    }))
    X_bg <- X_eval[bg_idx, , drop = FALSE]

    # drop zero-variance columns (avoid degeneracy)
    v2   <- matrixStats::colVars(as.matrix(X_eval))
    keep <- v2 > 0 & is.finite(v2)
    X_eval <- X_eval[, keep, drop = FALSE]
    X_bg   <- X_bg[,  keep, drop = FALSE]
    genes  <- colnames(X_eval)

    # cache path
    if (is.null(cache_dir)) cache_dir <- tools::R_user_dir("RFGeneRank", "cache")
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    key <- paste(
      digest::digest(serialize(fit, NULL), "xxhash64"),
      digest::digest(genes, "xxhash64"),
      digest::digest(X_bg, "xxhash64"),
      nsim, sep = "_"
    )
    cache_path <- file.path(cache_dir, paste0("shap_", key, ".rds"))

    # ---------- Robust predictor (auto-build a temp model if needed) ----------
    if (!is.null(pred_fun)) {
      pf_local <- match.fun(pred_fun)
      pred_wrapper <- function(object, newdata) pf_local(newdata)
    } else {
      # pull finalized learner + training features; if missing, build ephemeral model now
      mdl <- try(mdl0, silent = TRUE)
if (inherits(mdl, "try-error") || is.null(mdl)) {
  if (is.list(fit) && !is.null(fit$final_model)) mdl <- fit$final_model
}

feats <- try(fts0, silent = TRUE)
if (inherits(feats, "try-error") || is.null(feats)) {
  if (is.list(fit) && !is.null(fit$features)) feats <- fit$features
}
      if (is.null(mdl)) {
        # ---- build a temporary full-data model on the selected genes (not saved back) ----
        df_tr  <- as.data.frame(X_eval, check.names = FALSE)  # samples x genes
        trees <- try(pars$trees, silent = TRUE)
        if (inherits(trees, "try-error") || is.null(trees)) trees <- 1000L

        imp <- try(pars$importance, silent = TRUE)
        if (inherits(imp, "try-error") || is.null(imp)) imp <- "permutation"

        cw <- try(pars$class_weights, silent = TRUE)
        if (inherits(cw, "try-error")) cw <- NULL

        seed0 <- try(pars$seed, silent = TRUE)
        if (inherits(seed0, "try-error") || is.null(seed0)) seed0 <- seed
        mdl <- ranger::ranger(
          x = df_tr, y = y,
          num.trees     = trees,
          classification = TRUE,
          probability    = TRUE,
          importance     = imp,
          class.weights  = cw,
          seed           = seed0
        )
        feats <- colnames(df_tr)
      }
      if (!length(feats)) {
        stop("factor_dependence: training features are empty. The model must be trained on some genes.")
      }

      # map "COVID-19" -> "COVID.19" if needed
     pick_pos_col <- function(P, want) {
  cn <- colnames(P)
  if (is.null(cn)) stop("Model returned probs without column names.")
  if (want %in% cn) return(want)

  mn <- make.names(cn)
  if (want %in% mn) return(cn[match(want, mn)])

  details <- paste(cn, collapse = ", ")
  msg <- paste0("Positive label '", want, "' not found in predicted columns: ", details)
  stop(msg)
}

poslab <- pos

pred_wrapper <- function(object, newdata) {
  D <- as.data.frame(newdata)
  miss <- setdiff(feats, colnames(D))

  if (length(miss)) {
    details <- paste(utils::head(miss, 8), collapse = ", ") # include only a short preview of missing entries in the error details
    msg <- paste0(
      "newdata is missing ",
      length(miss),
      " trained features. Examples: ",
      details
    )
    stop(msg)
  }

  D <- D[, feats, drop = FALSE]

  if (inherits(mdl, "ranger")) {
    pr <- predict(mdl, data = D)
    P  <- pr$predictions
    if (is.null(P)) stop("Underlying ranger not probability-enabled; set probability=TRUE.")
    poscol <- pick_pos_col(P, poslab)
    return(as.numeric(P[, poscol, drop = TRUE]))
  }

  if (inherits(mdl, "randomForest")) {
    P <- predict(mdl, D, type = "prob")
    poscol <- pick_pos_col(P, poslab)
    return(as.numeric(P[, poscol, drop = TRUE]))
  }

  if (inherits(mdl, "train")) { # caret
    P <- predict(mdl, D, type = "prob")
    P <- as.data.frame(P)
    if (!(poslab %in% colnames(P))) {
      mn <- make.names(colnames(P))
      if (poslab %in% mn) poslab <- colnames(P)[match(poslab, mn)]
    }
    return(as.numeric(P[[poslab]]))
  }

  if (inherits(mdl, "glm")) {   # binomial
    P <- predict(mdl, newdata = D, type = "response")
    return(as.numeric(P))
  }

  details <- paste(class(mdl), collapse = ", ")
  msg <- paste0("Unsupported final_model class: ", details)
  stop(msg)
}
}
    # ---------------------------------------------------------------------------

    # quick sanity check
    p_bg <- try(pred_wrapper(fit, X_bg), silent = TRUE)
    if (inherits(p_bg, "try-error")) stop(conditionMessage(attr(p_bg, "condition")))
    if (sd(p_bg, na.rm = TRUE) < 1e-8) {
      stop("Background predictions are (near) constant; cannot compute meaningful SHAP. ",
           "Check that the model is probability-enabled and features match training.")
    }

    # cache or compute SHAP
    if (file.exists(cache_path)) {
      contr <- readRDS(cache_path)
    } else {
      contr <- fastshap::explain(
        object       = fit,
        X            = as.data.frame(X_bg),
        newdata      = as.data.frame(X_eval),
        nsim         = nsim,
        pred_wrapper = pred_wrapper
      )
      contr <- as.matrix(contr)
      saveRDS(contr, cache_path)
    }
  }

  ## ---------- Proxy path (fallback) ----------
  if (is.null(contr)) {
    used <- "proxy"
    # z-score per gene (genes x samples -> z, then transpose to samples x genes via 'contr' step)
    z <- t(scale(t(expr))); z[is.na(z)] <- 0

    imp <- imp0$importance
    names(imp) <- imp0$gene

    dirv <- if ("signed_importance" %in% names(imp0)) {
    d <- ifelse(imp0$direction > 0, 1, -1)
    names(d) <- imp0$gene
    d
    } else {
    setNames(rep(1, nrow(imp0)), imp0$gene)
    }

    g <- intersect(rownames(expr), names(imp))
    z <- z[g, , drop = FALSE]; imp <- imp[g]; dirv <- dirv[g]

    # samples x genes (IMPORTANT: correct orientation)
    contr <- t(z) * (imp * dirv)

    # keep only selected genes, preserving order
    common <- intersect(colnames(contr), genes)
    contr  <- contr[, common, drop = FALSE]
    X_eval <- X_eval[, common, drop = FALSE]
    genes  <- common
  }

  contr <- contr[, genes, drop = FALSE]   # ensure column order

  ## ---------- Per (gene x covariate) tests ----------
  res <- vector("list", length = length(genes) * length(covariates))
  k <- 0L

  for (cv in covariates) {
    if (!cv %in% colnames(cd)) next
    xcv <- cd[[cv]]

    for (gn in genes) {
      k <- k + 1L
      v <- contr[, gn]

      if (is.numeric(xcv)) {
        r_v <- resid(lm(v ~ y))
        r_x <- resid(lm(xcv ~ y))
        ct  <- stats::cor.test(r_v, r_x, method = "spearman", exact = FALSE)
        res[[k]] <- data.frame(
          gene = gn, covariate = cv, test = "partial Spearman",
          stat = unname(ct$estimate), pval = ct$p.value,
          effect = sign(unname(ct$estimate)),
          stringsAsFactors = FALSE
        )
      } else {
        xcvf <- droplevels(as.factor(xcv))
        fit_lm <- lm(v ~ y + xcvf)
        a <- anova(fit_lm)
        p <- tryCatch(a["xcvf", "Pr(>F)"], error = function(e) NA_real_)
        means <- tapply(v, xcvf, mean, na.rm = TRUE)
        rng <- if (all(is.finite(means))) diff(range(means)) else NA_real_
        res[[k]] <- data.frame(
          gene = gn, covariate = cv, test = "LM (v ~ y + covariate)",
          stat = rng, pval = p, effect = ifelse(is.na(rng), NA, sign(rng)),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, res)
  out$fdr <- p.adjust(out$pval, method = fdr_method)
  out$dependent <- out$fdr <= 0.05
  rownames(out) <- NULL
  attr(out, "method")     <- used
  attr(out, "assay")      <- "expr"
  attr(out, "genes_used") <- genes
  attr(out, "cache_file") <- cache_path
  out
}
