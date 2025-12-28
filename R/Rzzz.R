# ---- R/Rzzz.R ---------------------------------------------------------------

# This file is intentionally minimal - package hooks or onLoad can go here.
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("RFGeneRank: Cross-validated, stable predictive gene ranking for transcriptomics")
}

