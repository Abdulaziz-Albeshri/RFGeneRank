#' Accessors for GeneRankFit
#'
#' Getter/setter accessors for `GeneRankFit` slots.
#'
#' @param x A `GeneRankFit` object.
#' @param value Replacement value.
#'
#' @return
#' For getters:
#' \describe{
#'   \item{`params(x)`}{A list of training parameters and metadata.}
#'   \item{`oof(x)`}{A list of out-of-fold results including predictions and labels.}
#'   \item{`imp(x)`}{A `data.frame` containing aggregated variable importance scores.}
#'   \item{`calibration(x)`}{A list describing calibration settings.}
#' }
#'
#' For setters:
#' \describe{
#'   \item{`imp(x) <- value`}{Returns the updated `GeneRankFit` object with modified importance data.}
#'   \item{`calibration(x) <- value`}{Returns the updated `GeneRankFit` object with modified calibration settings.}
#' }
#'
#' @examples
#' fit <- methods::new("GeneRankFit")
#'
#' # getters
#' params(fit)
#' oof(fit)
#' imp(fit)
#' calibration(fit)
#'
#' # setters
#' imp(fit) <- data.frame(
#'   gene = c("GeneA", "GeneB"),
#'   importance = c(0.8, 0.3)
#' )
#' calibration(fit) <- list(
#'   method = "platt",
#'   fun = function(p) p
#' )
#'
#' imp(fit)
#' calibration(fit)
#'
#' @name GeneRankFit-accessors
NULL

# params
#' @rdname GeneRankFit-accessors
#' @export
setGeneric("params", function(x) standardGeneric("params"))
#' @rdname GeneRankFit-accessors
#' @export
setMethod("params", "GeneRankFit", function(x) .rfgr_params(x))

# oof
#' @rdname GeneRankFit-accessors
#' @export
setGeneric("oof", function(x) standardGeneric("oof"))
#' @rdname GeneRankFit-accessors
#' @export
setMethod("oof", "GeneRankFit", function(x) .rfgr_oof(x))

# imp
#' @rdname GeneRankFit-accessors
#' @export
setGeneric("imp", function(x) standardGeneric("imp"))
#' @rdname GeneRankFit-accessors
#' @export
setMethod("imp", "GeneRankFit", function(x) .rfgr_imp(x))


# imp (setter)
#' @rdname GeneRankFit-accessors
#' @export
setGeneric("imp<-", function(x, value) standardGeneric("imp<-"))
#' @rdname GeneRankFit-accessors
#' @export
setReplaceMethod("imp", "GeneRankFit", function(x, value) {
  .rfgr_imp_set(x, value)
})

# calibration (getter)
#' @rdname GeneRankFit-accessors
#' @export
setGeneric("calibration", function(x) standardGeneric("calibration"))
#' @rdname GeneRankFit-accessors
#' @export
setMethod("calibration", "GeneRankFit", function(x) .rfgr_calibration(x))

# calibration (setter)
#' @rdname GeneRankFit-accessors
#' @export
setGeneric("calibration<-", function(x, value) standardGeneric("calibration<-"))
#' @rdname GeneRankFit-accessors
#' @export
setReplaceMethod("calibration", "GeneRankFit", function(x, value) {
  .rfgr_calibration_set(x, value)
})