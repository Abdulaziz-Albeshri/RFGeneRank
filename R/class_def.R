# ---- R/class_def.R -----------------------------------------------------------
# S4 class that carries CV params, OOF predictions, importances,
# finalized model, and optional calibration function.

setClass(
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
    if (length(object@oof)) {
      if (is.null(object@oof$prob) || is.null(object@oof$y)) {
        msgs <- c(msgs, "'oof' must contain at least $prob and $y")
      } else {
        if (!is.matrix(object@oof$prob)) msgs <- c(msgs, "oof$prob must be a matrix")
        if (!is.factor(object@oof$y))   msgs <- c(msgs, "oof$y must be a factor")
      }
    }

    # calibration must be a list (possibly empty) with optional $method and $fun
    if (!is.list(object@calibration)) {
      msgs <- c(msgs, "'calibration' must be a list (use list() when unset)")
    } else if (!is.null(object@calibration$fun) && !is.function(object@calibration$fun)) {
      msgs <- c(msgs, "'calibration$fun' must be a function")
    }

    if (length(msgs)) msgs else TRUE
  }
)

# Pretty show method (minimal)
setMethod("show", "GeneRankFit", function(object) {
  cat("GeneRankFit\n")
  if (length(object@params)) {
    p <- object@params
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
  if (nrow(object@imp)) {
    cat(sprintf("  features ranked: %d\n", nrow(object@imp)))
  }
  if (!is.null(object@final_model)) {
    cat("  finalized model: PRESENT\n")
  } else {
    cat("  finalized model: <none>\n")
  }
  if (length(object@calibration) && is.character(object@calibration$method)) {
    cat(sprintf("  calibration: %s\n", object@calibration$method))
  } else {
    cat("  calibration: <none>\n")
  }
})

# small infix helper for show()
`%||%` <- function(a, b) if (is.null(a)) b else a
