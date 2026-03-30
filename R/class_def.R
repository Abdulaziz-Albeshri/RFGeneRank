# ---- R/class_def.R -----------------------------------------------------------
# S4 class that carries CV params, OOF predictions, importances,
# finalized model, and optional calibration function.

methods::setClass(
  "GeneRankFit",
  slots = c(
    params       = "list",         # training params & metadata
    oof          = "list",         # out-of-fold: prob (matrix), pred (factor), y (factor)
    imp          = "data.frame",   # aggregated importances
    features     = "character",    # ranked feature vector (gene IDs)
    final_model  = "ANY",          # ranger model after finalize()
    var_pca      = "numeric",      # variance explained (for info/plots)
    calibration  = "list"          # list(method, fun); use empty list() when not set
  ),
  prototype = list(
    params       = list(),
    oof          = list(),
    imp          = data.frame(),
    features     = character(),
    final_model  = NULL,
    var_pca      = numeric(),
    calibration  = list()          # <-- IMPORTANT: empty list, not NULL
  ),
  validity = function(object) {
  msgs <- character()

  # oof minimal structure check (lenient to allow partial fits)
  oof0 <- .rfgr_oof(object)
  if (length(oof0)) {
    if (is.null(oof0$prob) || is.null(oof0$y)) {
      msgs <- c(msgs, "'oof' must contain at least $prob and $y")
    } else {
      if (!is.matrix(oof0$prob)) msgs <- c(msgs, "oof$prob must be a matrix")
      if (!is.factor(oof0$y))   msgs <- c(msgs, "oof$y must be a factor")
    }
  }

  # calibration must be a list (possibly empty) with optional $method and $fun
  cal0 <- .rfgr_calibration(object)
  if (!is.list(cal0)) {
    msgs <- c(msgs, "'calibration' must be a list (use list() when unset)")
  } else if (!is.null(cal0$fun) && !is.function(cal0$fun)) {
    msgs <- c(msgs, "'calibration$fun' must be a function")
  }

  if (length(msgs)) msgs else TRUE
  }
)

# Pretty show method (minimal)
methods::setMethod("show", "GeneRankFit", function(object) {
  cat("GeneRankFit\n")

  pars <- .rfgr_params(object)
  if (length(pars)) {
    p <- pars
    cat(sprintf("  k=%s, trees=%s, importance=%s\n",
                as.character(p$k), as.character(p$trees), as.character(p$importance)))

    if (!is.null(p$fold_batch_correction) && isTRUE(p$fold_batch_correction)) {
      cat("  fold-safe batch correction: ON\n")
    }
    if (!is.null(p$filter_low_expr) && isTRUE(p$filter_low_expr)) {
      cat(sprintf("  CV-safe filtering: keep >= %.0f%% nonzero expr\n", 100*as.numeric(p$min_prop)))
    }
    if (!is.null(p$transform) && p$transform != "none") {
      cat(sprintf("  transform: %s\n", p$transform))
    }
    if (!is.null(p$standardize) && isTRUE(p$standardize)) {
      cat("  standardize: z-score by train mean/SD\n")
    }
  }

  imp0 <- .rfgr_imp(object)
  if (nrow(imp0)) {
    cat(sprintf("  features ranked: %d\n", nrow(imp0)))
  }

  mdl0 <- .rfgr_final_model(object)
  if (!is.null(mdl0)) {
    cat("  finalized model: PRESENT\n")
  } else {
    cat("  finalized model: <none>\n")
  }

  cal0 <- .rfgr_calibration(object)
  if (length(cal0) && is.character(cal0$method)) {
    cat(sprintf("  calibration: %s\n", cal0$method))
  } else {
    cat("  calibration: <none>\n")
  }
})


