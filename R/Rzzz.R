# ---- R/Rzzz.R ---------------------------------------------------------------

# This file is intentionally minimal - package hooks or onLoad can go here.
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("RFGeneRank: CV-stable predictive ranking for transcriptomics")
}
