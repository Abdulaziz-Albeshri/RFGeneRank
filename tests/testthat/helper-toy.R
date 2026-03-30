toy_se <- function(ngenes = 30, nsamp = 12, seed = 1) {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(seed)

  expr <- matrix(rnorm(ngenes * nsamp), nrow = ngenes)
  rownames(expr) <- paste0("g", seq_len(ngenes)) 
  colnames(expr) <- paste0("s", seq_len(nsamp))

  y <- factor(rep(c("Control", "Case"), each = nsamp/2))

  SummarizedExperiment::SummarizedExperiment(
    assays  = list(expr = expr),
    colData = data.frame(
      label = y,
      batch = rep(1:2, length.out = nsamp),
      row.names = colnames(expr)
    )
  )
}