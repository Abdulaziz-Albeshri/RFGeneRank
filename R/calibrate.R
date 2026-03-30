# ---- R/calibrate.R -----------------------------------------------------------
#' Calibrate out-of-fold probabilities (isotonic or Platt)
#'
#' Fits a calibration model on the out-of-fold (OOF) positive-class probabilities
#' and stores it inside the GeneRankFit object. No data leakage: uses OOF only.
#' @importFrom stats predict
#' @param fit GeneRankFit containing out-of-fold predictions and labels.
#' @param method "isotonic" or "platt"
#' @return GeneRankFit with stored calibration information.
#' @examples
#'
#' # Example probabilities
#' p <- runif(5)
#'
#' # Case 1: no calibration stored → returns original p
#' fit0 <- methods::new("GeneRankFit")
#' apply_calibration(fit0, p)
#'
#' # Case 2: simple calibration function (illustration)
#' fit1 <- methods::new("GeneRankFit")
#' calibration(fit1) <- list(fun = function(x) x^0.8)
#' apply_calibration(fit1, p)
#' @export
calibrate_oof <- function(fit, method = c("isotonic","platt")) {
  method <- match.arg(method)
  if (!inherits(fit, "GeneRankFit")) stop("fit must be GeneRankFit.")
  oof <- .rfgr_oof(fit)
  prob <- oof$prob
  y <- oof$y
  stopifnot(!is.null(prob), !is.null(y))
  sids <- rownames(prob)
  if (is.null(sids)) stop("oof$prob must have rownames (sample IDs).")
  if (is.null(names(y))) names(y) <- sids
  pos <- levels(y)[2]
  p   <- as.numeric(prob[sids, pos])
  lbl <- as.numeric(y[sids] == pos)
  ok  <- is.finite(p) & !is.na(lbl)
  p   <- p[ok]; lbl <- lbl[ok]

  if (length(unique(lbl)) < 2L) stop("Calibration needs both classes present in OOF.", call. = FALSE)

  if (method == "platt") {
    # logistic regression on scores
    df  <- data.frame(y = lbl, p = p)
    mdl <- stats::glm(y ~ p, family = stats::binomial(), data = df)
    cal_fun <- function(x) as.numeric(predict(mdl, newdata = data.frame(p = as.numeric(x)), type = "response"))
    fit <- .rfgr_calibration_set(fit, list(method = "platt", fun = cal_fun))
    return(fit)
  }

  # ---- Isotonic: avoid predict(isoreg) -- build a right-continuous step function ----
  ord <- order(p)
  p_s <- p[ord]
  y_s <- lbl[ord]
  iso <- stats::isoreg(p_s, y_s)   # monotone fit on sorted scores

  # iso$yf are fitted values at p_s; create a right-continuous step function
  step_fun <- stats::approxfun(
    x = p_s,
    y = iso$yf,
    method = "constant",
    f = 1,                   # right-continuous
    yleft  = min(iso$yf),
    yright = max(iso$yf),
    ties = "ordered"
  )

  cal_fun <- function(x) {
    z <- step_fun(as.numeric(x))
    z <- pmin(pmax(z, 0), 1)
    as.numeric(z)
  }

  fit <- .rfgr_calibration_set(fit, list(method = "isotonic", fun = cal_fun))
  fit
}

#' Apply stored calibration to a numeric vector of probabilities
#'
#' @param fit GeneRankFit with stored calibration information.
#' @param p numeric vector of positive-class probabilities
#' @return numeric vector of calibrated probabilities (or original if none stored)
#' @export
#' @examples
#'
#' # Example probabilities
#' p <- runif(5)
#'
#' # Minimal GeneRankFit without calibration (returns original p)
#' fit0 <- methods::new("GeneRankFit")
#' apply_calibration(fit0, p)
apply_calibration <- function(fit, p) {
  if (!inherits(fit, "GeneRankFit")) stop("fit must be GeneRankFit.")
  cal <- .rfgr_calibration(fit)
  if (is.list(cal) && is.function(cal$fun)) {
  return(cal$fun(p))
}
  p
}
