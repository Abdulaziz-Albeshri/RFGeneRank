# ---- R/compute_class_weights.R ----------------------------------------------
#' Compute normalized inverse-frequency class weights (internal)
#'
#' @description
#' Compute a named numeric vector of class weights where rarer classes
#' receive larger weights. Weights are normalized to have mean 1.
#'
#' Examples:
#' y <- c("A","A","A","B","B")
#' .compute_class_weights(y)
#'
#' @keywords internal
#' @noRd
.compute_class_weights <- function(y, na.rm = TRUE) {
  if (is.null(y)) stop("`.compute_class_weights()`: `y` is NULL.")
  if (na.rm) y <- y[!is.na(y)]
  if (!is.factor(y)) y <- factor(y)

  tab <- table(y)
  if (any(tab == 0)) {
    # Shouldn't happen for a factor built from y, but be safe
    tab <- tab[tab > 0]
  }
  if (length(tab) < 2L) {
    stop("`.compute_class_weights()` needs at least two classes.")
  }

  # Inverse-frequency relative to the median class count
  w <- as.numeric(stats::median(tab) / tab)
  names(w) <- names(tab)

  # Normalize so mean weight = 1 (stable across class counts)
  w / mean(w)
}
