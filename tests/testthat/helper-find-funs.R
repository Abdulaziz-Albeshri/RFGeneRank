
.rfgr_ns <- function() asNamespace("RFGeneRank")

.rfgr_get_ns <- function(name) {
  ns <- .rfgr_ns()
  if (!exists(name, envir = ns, inherits = FALSE)) return(NULL)
  get(name, envir = ns, inherits = FALSE)
}

.rfgr_get_any <- function(candidates) {
  for (nm in candidates) {
    fn <- .rfgr_get_ns(nm)
    if (!is.null(fn)) return(fn)
  }
  NULL
}

.rfgr_skip_if_missing <- function(x, msg) {
  if (is.null(x)) testthat::skip(msg)
}