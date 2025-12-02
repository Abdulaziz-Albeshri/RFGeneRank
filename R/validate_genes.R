# ---- R/validate_genes.R ----------------------------------------------------
#' Validate a ranked gene set with alternate learners via k-fold CV
#'
#' @param se A SummarizedExperiment whose assay is genes-by-samples.
#' @param genes Character vector of gene IDs (must match \code{rownames(assay(se))}).
#' @param methods Character vector; any of \code{c("ranger","glmnet","xgboost")}.
#' @param k Integer number of folds (default 5).
#' @param seed Integer RNG seed for reproducibility.
#' @param model_params Named list of per-method parameter lists,
#'   e.g. \code{list(
#'     ranger  = list(num.trees=1000, importance="none", nthread=1),
#'     glmnet  = list(alpha=0.5, standardize="zscore", clip=8, eps_sd=1e-8,
#'                    use_class_weights=TRUE, lambda="lambda.min"),
#'     xgboost = list(nrounds=400, eta=0.05, max_depth=4, nthread=1)
#'   )}.
#' @param label_col Column name in \code{colData(se)} with the class labels.
#'   Must be a binary factor (exactly \code{2} levels).
#' @param positive Optional; the positive class label. If \code{NULL}, uses \code{levels(y)[2]}.
#' @param calibrate One of \code{c("platt","isotonic","none")}; applied per fold.
#' @param thr_metric One of \code{c("youden","f1","cost")} for threshold selection.
#' @param cost Named numeric vector \code{c(fp=1, fn=1)} if \code{thr_metric == "cost"}.
#' @return A list with \code{$summary} and one entry per method. Each entry contains:
#'   \code{$oof} (data.frame with p_raw, p_cal, p_use, y, fold),
#'   \code{$calibrator_requested}, \code{$calibrator_used}, \code{$threshold},
#'   \code{$conf_mat}, and \code{$metrics} with both diagnostics (raw/cal)
#'   and final guardrailed metrics (\code{auc_final}, \code{ece_final}, \code{brier_final}).
#' @examples
#' # Toy expression matrix: 10 genes x 12 samples
#' expr <- matrix(
#'   stats::rnorm(10 * 12),
#'   nrow = 10,
#'   dimnames = list(
#'     paste0("gene", 1:10),
#'     paste0("sample", 1:12)
#'   )
#' )
#'
#' # Binary labels as a factor (balanced)
#' label <- factor(rep(c("A", "B"), each = 6))
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays  = list(expr = expr),
#'   colData = S4Vectors::DataFrame(label = label)
#' )
#'
#' # Use the first 5 genes as a toy signature
#' genes <- rownames(expr)[1:5]
#'
#' # Validate using ranger only (fast and robust for examples)
#' res <- validate_genes(
#'   se        = se,
#'   genes     = genes,
#'   methods   = "ranger",
#'   k         = 3,
#'   seed      = 1,
#'   label_col = "label"
#' )
#'
#' # Inspect out-of-fold AUROC
#' res$summary$auc_final
#' @export
validate_genes <- function(
  se, genes,
  methods = c("glmnet","xgboost","ranger"),
  k = 5, seed = 1L,
  model_params = list(),
  label_col = "state",
  positive = NULL,
  calibrate = c("platt","isotonic","none"),
  thr_metric = c("youden","f1","cost"),
  cost = c(fp = 1, fn = 1)
) {
  # ---- input checks
  stopifnot(inherits(se, "SummarizedExperiment"))
  calibrate  <- match.arg(calibrate)
  thr_metric <- match.arg(thr_metric)
  methods    <- unique(methods)

  Xall <- SummarizedExperiment::assay(se)
  if (missing(genes) || is.null(genes) || !length(genes)) {
    stop("Please provide a non-empty character vector `genes` (e.g., top_genes(fit, n=100)$gene).")
  }
  if (!all(genes %in% rownames(Xall))) {
    miss <- setdiff(genes, rownames(Xall))
    stop("Missing genes in assay(se): ", paste(miss, collapse = ", "))
  }
  X <- Xall[genes, , drop = FALSE]

  y <- SummarizedExperiment::colData(se)[[label_col]]
  if (!is.factor(y)) y <- factor(y)
  y <- droplevels(y)
  if (nlevels(y) != 2L) stop("`label_col` must be a binary factor.")
  if (is.null(positive)) positive <- levels(y)[2]
  if (!(positive %in% levels(y))) stop("`positive` not found in levels(y).")
  # Ensure second level is the positive class for consistent AUC orientation
  if (levels(y)[2] != positive) {
    y <- stats::relevel(y, ref = setdiff(levels(y), positive))
  }
  lev <- levels(y)  # lev[2] is positive

  
  folds <- caret::createFolds(y, k = k, returnTrain = TRUE)
  if (!length(folds)) stop("Could not create CV folds.")

  # ---- helpers --------------------------------------------------------------
  `%||%` <- function(x, y) if (is.null(x)) y else x

  compute_class_weights <- function(y) {
    p <- table(y)/length(y)
    w <- 1/as.numeric(p)
    stats::setNames(as.numeric(w/mean(w)), names(p))
  }

  # Expected Calibration Error (binary)
  ece_metric <- function(p, y, bins = 10) {
    z  <- as.integer(y == levels(y)[2])
    br <- cut(p, breaks = seq(0, 1, length.out = bins + 1), include.lowest = TRUE)
    by <- split(data.frame(p, z), br)
    w  <- vapply(by, nrow, 0.0)
    if (sum(w) == 0) return(NA_real_)
    w  <- w / sum(w)
    conf <- vapply(by, function(df) if (nrow(df)) mean(df$p) else 0, 0.0)
    acc  <- vapply(by, function(df) if (nrow(df)) mean(df$z) else 0, 0.0)
    sum(w * abs(conf - acc))
  }

  # Brier score (binary)
  brier <- function(p, y) {
    z <- as.integer(y == levels(y)[2])
    mean((p - z)^2)
  }

  # Isotonic calibration (fold-wise), warning-free
  isotonic_calibrator <- function(p_tr, y_tr) {
    z   <- as.integer(y_tr == levels(y_tr)[2])
    iso <- stats::isoreg(p_tr, z)
    ord <- order(iso$x)
    x   <- iso$x[ord]
    yv  <- iso$yf[ord]
    keep <- !duplicated(x)
    ffun <- stats::approxfun(x[keep], yv[keep],
                             method = "constant", f = 1, rule = 2)
    function(p) {
      out <- ffun(p)
      pmin(1, pmax(0, out))
    }
  }

  # Robust Platt calibration with auto-fallback to isotonic
  platt_calibrator <- function(p_tr, y_tr) {
    z <- as.integer(y_tr == lev[2])
    eps <- 1e-6
    ltr <- stats::qlogis(pmin(pmax(p_tr, eps), 1 - eps))

    warn_flag <- FALSE
    m <- try(withCallingHandlers({
      stats::glm(z ~ ltr, family = stats::binomial(),
                 control = stats::glm.control(maxit = 200))
    }, warning = function(w) {
      warn_flag <- TRUE
      invokeRestart("muffleWarning")
    }), silent = TRUE)

    if (!inherits(m, "try-error") && !warn_flag && all(is.finite(stats::coef(m)))) {
      function(p) {
        lp <- stats::qlogis(pmin(pmax(p, eps), 1 - eps))
        stats::plogis(stats::predict(m, newdata = data.frame(ltr = lp), type = "link"))
      }
    } else {
      iso <- isotonic_calibrator(p_tr, y_tr)
      # small identity blend prevents 0/1 collapse if knots are sparse
      function(p) {
        q <- iso(p)
        alpha <- 0.1
        pmin(1, pmax(0, (1 - alpha) * q + alpha * p))
      }
    }
  }

  # Threshold selection on provided probabilities
  choose_threshold <- function(p, y, metric = thr_metric, cost = cost) {
    roc <- pROC::roc(response = y, predictor = p, levels = rev(lev), quiet = TRUE)
    if (metric == "youden") {
      as.numeric(pROC::coords(roc, "best", best.method = "youden", ret = "threshold"))
    } else if (metric == "f1") {
      cand <- as.numeric(pROC::coords(roc, x = "threshold", ret = "threshold"))
      f1s <- vapply(cand, function(t) {
        levs <- levels(y)
        pr  <- factor(ifelse(p >= t, levs[2], levs[1]), levels = levs)
        tab <- table(True = y, Pred = pr)
        tp <- tab[levs[2], levs[2]]; fp <- tab[levs[1], levs[2]]
        fn <- tab[levs[2], levs[1]]
        prec <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
        rec  <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
        if (prec + rec == 0) 0 else 2 * prec * rec / (prec + rec)
      }, 0.0)
      cand[which.max(f1s)]
    } else { # cost-based
      cand <- seq(0.01, 0.99, by = 0.01)
      losses <- vapply(cand, function(t) {
        levs <- levels(y)
        pr  <- factor(ifelse(p >= t, levs[2], levs[1]), levels = levs)
        tab <- table(True = y, Pred = pr)
        fp <- tab[levs[1], levs[2]]; fn <- tab[levs[2], levs[1]]
        cost["fp"] * fp + cost["fn"] * fn
      }, 0.0)
      cand[which.min(losses)]
    }
  }

  # ---------------- fit/predict per method (returns p_tr, p_te) ----------------
  fit_predict <- function(method, Xtr, ytr, Xte, params) {

    if (method == "glmnet") {
      # ---- fold-safe, configurable standardization ----
      std_style <- params$standardize %||% "zscore"   # "zscore","center","robust","none"
      clip_val  <- params$clip       %||% NULL
      eps_sd    <- params$eps_sd     %||% 1e-8

      # compute train stats and apply
      if (std_style == "zscore") {
        mu <- matrixStats::rowMeans2(Xtr)
        sd <- matrixStats::rowSds(Xtr); sd[sd < eps_sd] <- 1
        scale_apply <- function(M) t((t(M) - mu)/sd)
      } else if (std_style == "center") {
        mu <- matrixStats::rowMeans2(Xtr)
        scale_apply <- function(M) t(t(M) - mu)
      } else if (std_style == "robust") {
        med <- matrixStats::rowMedians(Xtr)
        mad <- matrixStats::rowMads(Xtr); mad[mad < eps_sd] <- 1
        scale_apply <- function(M) t((t(M) - med)/mad)
      } else {
        scale_apply <- function(M) M
      }
      Xtr_s <- scale_apply(Xtr)
      Xte_s <- scale_apply(Xte)

      # optional clipping to limit extreme z-scores
      if (!is.null(clip_val)) {
        Xtr_s[Xtr_s >  clip_val] <-  clip_val; Xtr_s[Xtr_s < -clip_val] <- -clip_val
        Xte_s[Xte_s >  clip_val] <-  clip_val; Xte_s[Xte_s < -clip_val] <- -clip_val
      }

      # class weights (optional, improves imbalance handling)
      ybin <- as.integer(ytr == lev[2])
      w_tr <- if (isTRUE(params$use_class_weights)) {
        cw <- table(ytr); cw <- as.numeric(mean(cw)/cw)
        ifelse(ytr == lev[1], cw[1], cw[2])
      } else NULL

      alpha <- params$alpha %||% 0.5
      fam   <- "binomial"
      lam_pick <- params$lambda %||% "lambda.min"  # or "lambda.1se"

      cv <- glmnet::cv.glmnet(
        x = t(Xtr_s), y = ybin,
        alpha = alpha, family = fam,
        type.measure = "auc",
        standardize = FALSE,    # Standardization is performed manually as part of the fold-safe procedure
        weights = w_tr
      )
      b  <- stats::coef(cv, s = lam_pick)

      lin_tr <- drop(Matrix::t(Xtr_s) %*% b[-1]) + b[1]
      p_tr   <- stats::plogis(lin_tr)
      lin_te <- drop(Matrix::t(Xte_s) %*% b[-1]) + b[1]
      p_te   <- stats::plogis(lin_te)
      return(list(p_tr = as.numeric(p_tr), p_te = as.numeric(p_te)))

    } else if (method == "xgboost") {
      # gbtree: no standardization needed
      ybin <- as.integer(ytr == lev[2])

      # (optional) instance weights for imbalance
      if (isTRUE(params$use_class_weights)) {
        cw <- table(ytr); w_pos <- as.numeric(mean(cw)/cw[lev[2]])
        w_neg <- as.numeric(mean(cw)/cw[lev[1]])
        w_tr <- ifelse(ytr == lev[2], w_pos, w_neg)
      } else w_tr <- NULL

      dtr  <- xgboost::xgb.DMatrix(data = t(Xtr), label = ybin, weight = w_tr)
      dte  <- xgboost::xgb.DMatrix(data = t(Xte))

      param <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta        = params$eta        %||% 0.05,
        max_depth  = params$max_depth  %||% 4,
        subsample  = params$subsample  %||% 1,
        colsample_bytree = params$colsample_bytree %||% 1,
        nthread    = params$nthread    %||% 1
      )
      bst <- xgboost::xgb.train(params = param, data = dtr,
                                nrounds = params$nrounds %||% 400, verbose = 0)
      p_tr <- stats::predict(bst, dtr)
      p_te <- stats::predict(bst, dte)
      return(list(p_tr = as.numeric(p_tr), p_te = as.numeric(p_te)))

    } else if (method == "ranger") {
      df_tr <- data.frame(t(Xtr)); df_tr$y <- ytr
      cw <- compute_class_weights(ytr)  # inverse-prevalence
      mdl <- ranger::ranger(
        y ~ ., data = df_tr,
        probability   = TRUE,
        num.trees     = params$num.trees     %||% 1000,
        mtry          = params$mtry          %||% max(1L, floor(sqrt(nrow(Xtr)))),
        min.node.size = params$min.node.size %||% 1L,
        importance    = params$importance    %||% "none",
        class.weights = cw,
        num.threads   = params$nthread       %||% 1,
        seed          = NULL
      )
      p_tr <- predict(mdl, data = data.frame(t(Xtr)))$predictions[, lev[2]]
      p_te <- predict(mdl, data = data.frame(t(Xte)))$predictions[, lev[2]]
      return(list(p_tr = as.numeric(p_tr), p_te = as.numeric(p_te)))

    } else {
      stop("Unsupported method: ", method)
    }
  }

  # ---- run one method across folds -----------------------------------------
  run_one_method <- function(meth) {
    params <- model_params[[meth]] %||% list()
    n <- length(y)
    p_raw <- numeric(n); p_cal <- numeric(n); fold_id <- integer(n)

    for (i in seq_along(folds)) {
      tr <- folds[[i]]; te <- setdiff(seq_along(y), tr)

      pr <- fit_predict(meth, X[, tr, drop = FALSE], y[tr],
                        X[, te, drop = FALSE], params)

      p_raw[te]   <- pr$p_te
      fold_id[te] <- i

      if (calibrate == "platt") {
        cal <- platt_calibrator(pr$p_tr, y[tr])
        p_cal[te] <- cal(pr$p_te)
      } else if (calibrate == "isotonic") {
        cal <- isotonic_calibrator(pr$p_tr, y[tr])
        p_cal[te] <- cal(pr$p_te)
      } else {
        p_cal[te] <- pr$p_te
      }
    }

    # ---- OOF diagnostics: raw vs calibrated
    roc_raw <- pROC::roc(response = y, predictor = p_raw, levels = rev(lev), quiet = TRUE)
    roc_cal <- pROC::roc(response = y, predictor = p_cal, levels = rev(lev), quiet = TRUE)

    auc_raw <- as.numeric(pROC::auc(roc_raw))
    auc_cal <- as.numeric(pROC::auc(roc_cal))
    ece_raw <- ece_metric(p_raw, y)
    ece_cal <- ece_metric(p_cal, y)
    br_raw  <- brier(p_raw, y)
    br_cal  <- brier(p_cal, y)

    # ---- Guardrail: auto-revert if calibration harms
    harm_auc <- (auc_cal < (auc_raw - 0.02))
    harm_ece <- (ece_cal > (ece_raw + 0.05))
    if (harm_auc || harm_ece) {
      p_use <- p_raw
      calibrator_used <- "none (auto-revert)"
    } else {
      p_use <- p_cal
      calibrator_used <- calibrate
    }

    # ---- Final metrics computed on the probabilities used in downstream evaluation
    roc_use <- pROC::roc(response = y, predictor = p_use, levels = rev(lev), quiet = TRUE)
    auc_use <- as.numeric(pROC::auc(roc_use))
    ece_use <- ece_metric(p_use, y)
    br_use  <- brier(p_use, y)

    # threshold on "used" probabilities
    thr <- choose_threshold(p_use, y, thr_metric, cost)

    pred <- factor(ifelse(p_use >= thr, lev[2], lev[1]), levels = lev)
    tab  <- table(True = y, Pred = pred)

    tp <- tab[lev[2], lev[2]]; tn <- tab[lev[1], lev[1]]
    fp <- tab[lev[1], lev[2]]; fn <- tab[lev[2], lev[1]]
    acc <- (tp + tn) / sum(tab)
    sens <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
    spec <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
    ppv  <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
    npv  <- if ((tn + fn) == 0) NA_real_ else tn / (tn + fn)

    out <- list(
      method = meth,
      oof = data.frame(
        p_raw = p_raw,
        p_cal = p_cal,
        p_use = p_use,
        y     = y,
        fold  = fold_id,
        stringsAsFactors = FALSE
      ),
      calibrator_requested = calibrate,
      calibrator_used      = calibrator_used,
      threshold  = thr,
      conf_mat   = tab,
      metrics = list(
        # diagnostics
        auc_raw = auc_raw, ece_raw = ece_raw, brier_raw = br_raw,
        auc_cal = auc_cal, ece_cal = ece_cal, brier_cal = br_cal,
        # final (guardrailed)
        auc_final = auc_use, ece_final = ece_use, brier_final = br_use,
        acc = acc, sens = sens, spec = spec, ppv = ppv, npv = npv
      )
    )
    class(out) <- "RFGeneRankFitLike"
    out
  }

  # run and summarize
  allowed <- c("ranger","glmnet","xgboost")
  unknown <- setdiff(methods, allowed)
  if (length(unknown)) stop("Unsupported methods: ", paste(unknown, collapse = ", "))

  fits <- lapply(methods, run_one_method)
  names(fits) <- methods

  summary <- data.frame(
    method       = methods,
    auc_final    = vapply(fits, function(z) z$metrics$auc_final,    0.0),
    ece_final    = vapply(fits, function(z) z$metrics$ece_final,    0.0),
    brier_final  = vapply(fits, function(z) z$metrics$brier_final,  0.0),
    threshold    = vapply(fits, function(z) z$threshold,            0.0),
    used_calib   = vapply(fits, function(z) z$calibrator_used,      character(1)),
    stringsAsFactors = FALSE
  )

  c(list(summary = summary), fits)
}
