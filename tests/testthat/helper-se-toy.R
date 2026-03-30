toy_se_cov <- function(ngenes = 120, nsamp = 20, seed = 1, entrez = FALSE) {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(seed)

  expr <- matrix(rnorm(ngenes * nsamp), nrow = ngenes)
  if (entrez) {
    rownames(expr) <- as.character(seq_len(ngenes))  # ENTRES-like: "1","2",...
  } else {
    rownames(expr) <- paste0("g", seq_len(ngenes))
  }
  colnames(expr) <- paste0("s", seq_len(nsamp))

  cd <- data.frame(
    state = factor(rep(c("Control","Case"), length.out = nsamp)),
    batch = factor(rep(c("B1","B2"), length.out = nsamp)),
    sex   = factor(rep(c("M","F"), length.out = nsamp)),
    row.names = colnames(expr)
  )

  SummarizedExperiment::SummarizedExperiment(
    assays  = list(expr = expr),
    colData = cd
  )
}