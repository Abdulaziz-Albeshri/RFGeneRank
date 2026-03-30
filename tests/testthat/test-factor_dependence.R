test_that("factor_dependence (proxy) runs and returns a result table", {
  skip_if_not_installed("SummarizedExperiment")

  set.seed(1)

  ngenes <- 120
  nsamp  <- 20

  expr <- matrix(rnorm(ngenes * nsamp), nrow = ngenes)
  rownames(expr) <- paste0("g", 1:ngenes)
  colnames(expr) <- paste0("s", 1:nsamp)

  cd <- data.frame(
    state = factor(rep(c("Control","Case"), each = nsamp/2)), 
    batch = factor(rep(c("B1","B2"), length.out = nsamp)),
    sex   = factor(rep(c("M","F"), length.out = nsamp)),
    row.names = colnames(expr)
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(expr = expr),
    colData = cd
  )

  fit <- RFGeneRank::rank_genes(
    se = se,
    label_col = "state",   
    batch_col = "batch",
    trees = 50,
    n_top = 50
  )

  res <- RFGeneRank::factor_dependence(
    fit = fit,
    se = se,
    covariates = "sex",
    method = "proxy",
    ngenes = 30,
    gene_selection = "importance",
    seed = 1
  )

  expect_true(is.data.frame(res) || is.list(res))

  if (is.data.frame(res)) {
    expect_true(nrow(res) >= 1)
  }
})

test_that("factor_dependence errors when assay expr is missing", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(rnorm(20), nrow = 5)),
    colData = S4Vectors::DataFrame(state = factor(c("A", "A", "B", "B")))
  )

  prob <- matrix(runif(4), ncol = 1, dimnames = list(colnames(SummarizedExperiment::assay(se)), "B"))
  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = factor(c("A", "A", "B", "B")))
  fit@imp <- data.frame(gene = paste0("g", 1:5), importance = 1:5)

  expect_error(
    RFGeneRank::factor_dependence(fit, se, covariates = "x", method = "proxy"),
    "Assay 'expr' not found"
  )
})

test_that("factor_dependence errors on OOF sample mismatch with expr columns", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  expr <- matrix(rnorm(6 * 4), nrow = 6)
  rownames(expr) <- paste0("g", 1:6)
  colnames(expr) <- paste0("s", 1:4)

  cd <- S4Vectors::DataFrame(state = factor(c("A", "A", "B", "B")))
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd)

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(
    prob = matrix(runif(4), ncol = 1, dimnames = list(c("s1", "s2", "s3", "BAD"), "B")),
    y = factor(c("A", "A", "B", "B"))
  )
  fit@imp <- data.frame(gene = rownames(expr), importance = 1:6)

  expect_error(
    RFGeneRank::factor_dependence(fit, se, covariates = "x", method = "proxy"),
    "Sample ID mismatch: some OOF sample IDs not in assay 'expr'"
  )
})

test_that("factor_dependence errors when state is missing or not binary", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  expr <- matrix(rnorm(6 * 6), nrow = 6)
  rownames(expr) <- paste0("g", 1:6)
  colnames(expr) <- paste0("s", 1:6)

  cd1 <- S4Vectors::DataFrame(sex = factor(rep(c("M", "F"), 3)))
  rownames(cd1) <- colnames(expr)
  se1 <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd1)

  fit1 <- methods::new("GeneRankFit")
  fit1@oof <- list(
    prob = matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B")),
    y = factor(rep(c("A", "B"), 3))
  )
  fit1@imp <- data.frame(gene = rownames(expr), importance = 1:6)

  expect_error(
    RFGeneRank::factor_dependence(fit1, se1, covariates = "sex", method = "proxy"),
    "colData\\(se\\)\\$state is required"
  )

  cd2 <- S4Vectors::DataFrame(state = factor(c("A", "B", "C", "A", "B", "C")))
  rownames(cd2) <- colnames(expr)
  se2 <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd2)

  fit2 <- methods::new("GeneRankFit")
  fit2@oof <- list(
    prob = matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B")),
    y = factor(c("A", "B", "C", "A", "B", "C"))
  )
  fit2@imp <- data.frame(gene = rownames(expr), importance = 1:6)

  expect_error(
    RFGeneRank::factor_dependence(fit2, se2, covariates = "state", method = "proxy"),
    "Binary classification expected"
  )
})

test_that("factor_dependence validates positive argument branches", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  expr <- matrix(rnorm(6 * 6), nrow = 6)
  rownames(expr) <- paste0("g", 1:6)
  colnames(expr) <- paste0("s", 1:6)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 3)),
    age = c(20, 30, 40, 50, 60, 70)
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd)

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(
    prob = matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B")),
    y = cd$state
  )
  fit@imp <- data.frame(gene = rownames(expr), importance = 1:6)

  expect_error(
    RFGeneRank::factor_dependence(
      fit, se, covariates = "age",
      method = "proxy", positive = 3
    ),
    "`positive` must be NULL"
  )

  expect_error(
    RFGeneRank::factor_dependence(
      fit, se, covariates = "age",
      method = "proxy", positive = "Z"
    ),
    "Chosen positive class"
  )
})

test_that("factor_dependence errors when importance table is unavailable for importance selection", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  expr <- matrix(rnorm(6 * 6), nrow = 6)
  rownames(expr) <- paste0("g", 1:6)
  colnames(expr) <- paste0("s", 1:6)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 3)),
    sex = factor(rep(c("M", "F"), 3))
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd)

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(
    prob = matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B")),
    y = cd$state
  )

  expect_error(
    RFGeneRank::factor_dependence(
      fit, se,
      covariates = "sex",
      method = "proxy",
      gene_selection = "importance"
    ),
    "requires an importance table"
  )
})

test_that("factor_dependence errors when variance selection finds no variable genes", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  expr <- matrix(1, nrow = 5, ncol = 4)
  rownames(expr) <- paste0("g", 1:5)
  colnames(expr) <- paste0("s", 1:4)

  cd <- S4Vectors::DataFrame(
    state = factor(c("A", "A", "B", "B")),
    age = c(20, 30, 40, 50)
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd)

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(
    prob = matrix(runif(4), ncol = 1, dimnames = list(colnames(expr), "B")),
    y = cd$state
  )
  fit@imp <- data.frame(gene = rownames(expr), importance = 1:5)

  expect_error(
    RFGeneRank::factor_dependence(
      fit, se,
      covariates = "age",
      method = "proxy",
      gene_selection = "variance",
      ngenes = 3
    ),
    "No variable genes found in assay 'expr'"
  )
})

test_that("factor_dependence proxy path returns tabular result with numeric and factor covariates", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  set.seed(1)

  expr <- matrix(rnorm(20 * 10), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:10)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 5)),
    sex = factor(rep(c("M", "F"), length.out = 10)),
    age = seq_len(10)
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd)

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(
    prob = matrix(runif(10), ncol = 1, dimnames = list(colnames(expr), "B")),
    y = cd$state
  )
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = seq_len(nrow(expr))
  )

  res <- RFGeneRank::factor_dependence(
    fit, se,
    covariates = c("sex", "age"),
    method = "proxy",
    ngenes = 5,
    gene_selection = "importance",
    positive = 2
  )

  expect_true(is.data.frame(res))
  expect_true(all(c("gene", "covariate", "test", "stat", "pval", "effect", "fdr", "dependent") %in% names(res)))
  expect_true(any(res$test == "LM (v ~ y + covariate)"))
  expect_true(any(res$test == "partial Spearman"))
  expect_identical(attr(res, "method"), "proxy")
  expect_identical(attr(res, "assay"), "expr")
  expect_true(is.character(attr(res, "genes_used")))
})

test_that("factor_dependence drops state from covariates but still works with remaining covariates", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  set.seed(2)

  expr <- matrix(rnorm(10 * 6), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:6)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 3)),
    sex = factor(c("M", "F", "M", "F", "M", "F"))
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd)

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(
    prob = matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B")),
    y = cd$state
  )
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = rev(seq_len(nrow(expr)))
  )

  res <- RFGeneRank::factor_dependence(
    fit, se,
    covariates = c("state", "sex"),
    method = "proxy",
    ngenes = 3,
    gene_selection = "importance"
  )

  expect_true(is.data.frame(res))
  expect_false("state" %in% unique(res$covariate))
  expect_true("sex" %in% unique(res$covariate))
})

test_that("factor_dependence with only state covariate returns a non-data-frame empty-like object", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  set.seed(3)

  expr <- matrix(rnorm(8 * 6), nrow = 8)
  rownames(expr) <- paste0("g", 1:8)
  colnames(expr) <- paste0("s", 1:6)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 3))
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd)

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(
    prob = matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B")),
    y = cd$state
  )
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = seq_len(nrow(expr))
  )

  res <- RFGeneRank::factor_dependence(
    fit, se,
    covariates = "state",
    method = "proxy",
    ngenes = 3,
    gene_selection = "importance"
  )

  expect_false(is.data.frame(res))
  expect_true(is.list(res) || is.null(res))
})

test_that("factor_dependence errors when top-importance genes do not overlap expr rownames", {
  expr <- matrix(rnorm(6 * 4), nrow = 6)
  rownames(expr) <- paste0("g", 1:6)
  colnames(expr) <- paste0("s", 1:4)

  cd <- S4Vectors::DataFrame(
    state = factor(c("A", "A", "B", "B")),
    sex   = factor(c("M", "F", "M", "F"))
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(4), ncol = 1, dimnames = list(colnames(expr), "B"))

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = cd$state)
  fit@imp <- data.frame(
    gene = c("x1", "x2", "x3"),
    importance = c(3, 2, 1)
  )

  expect_error(
    RFGeneRank::factor_dependence(
      fit, se,
      covariates = "sex",
      method = "proxy",
      ngenes = 2,
      gene_selection = "importance"
    ),
    "No overlap between top-importance genes and 'expr' rownames"
  )
})

test_that("factor_dependence errors when variance selection finds no variable genes", {
  expr <- matrix(1, nrow = 5, ncol = 4)
  rownames(expr) <- paste0("g", 1:5)
  colnames(expr) <- paste0("s", 1:4)

  cd <- S4Vectors::DataFrame(
    state = factor(c("A", "A", "B", "B")),
    age   = c(20, 30, 40, 50)
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(4), ncol = 1, dimnames = list(colnames(expr), "B"))

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = cd$state)

  expect_error(
    RFGeneRank::factor_dependence(
      fit, se,
      covariates = "age",
      method = "proxy",
      ngenes = 3,
      gene_selection = "variance"
    ),
    "No variable genes found in assay 'expr'"
  )
})

test_that("factor_dependence proxy path uses unsigned importance when signed_importance is absent", {
  set.seed(1)

  expr <- matrix(rnorm(8 * 6), nrow = 8)
  rownames(expr) <- paste0("g", 1:8)
  colnames(expr) <- paste0("s", 1:6)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 3)),
    sex   = factor(c("M", "F", "M", "F", "M", "F")),
    age   = c(25, 31, 28, 45, 38, 50)
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B"))

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = cd$state)
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = seq_len(nrow(expr))
  )

  res <- RFGeneRank::factor_dependence(
    fit, se,
    covariates = c("sex", "age"),
    method = "proxy",
    ngenes = 4,
    gene_selection = "importance",
    positive = 2
  )

  expect_true(is.data.frame(res))
  expect_true(all(c("gene", "covariate", "test", "stat", "pval", "effect", "fdr", "dependent") %in% names(res)))
  expect_identical(attr(res, "method"), "proxy")
  expect_identical(attr(res, "assay"), "expr")
  expect_true(is.character(attr(res, "genes_used")))
  expect_true(length(attr(res, "genes_used")) == 4)
  expect_true(is.character(attr(res, "cache_file")) || is.na(attr(res, "cache_file")))
  expect_true(any(res$test == "LM (v ~ y + covariate)"))
  expect_true(any(res$test == "partial Spearman"))
})

test_that("factor_dependence proxy path accepts positive = 1 and drops state from covariates", {
  set.seed(2)

  expr <- matrix(rnorm(10 * 6), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:6)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 3)),
    sex   = factor(c("M", "F", "M", "F", "M", "F"))
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "A"))

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = cd$state)
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = rev(seq_len(nrow(expr)))
  )

  res <- RFGeneRank::factor_dependence(
    fit, se,
    covariates = c("state", "sex"),
    method = "proxy",
    ngenes = 5,
    gene_selection = "importance",
    positive = 1
  )

  expect_true(is.data.frame(res))
  expect_false("state" %in% unique(res$covariate))
  expect_true("sex" %in% unique(res$covariate))
})

test_that("factor_dependence proxy path uses signed_importance branch when available", {
  set.seed(3)

  expr <- matrix(rnorm(6 * 6), nrow = 6)
  rownames(expr) <- paste0("g", 1:6)
  colnames(expr) <- paste0("s", 1:6)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 3)),
    sex   = factor(rep(c("M", "F"), 3))
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B"))

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = cd$state)
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = c(6, 5, 4, 3, 2, 1),
    signed_importance = c(6, -5, 4, -3, 2, -1),
    direction = c(1, -1, 1, -1, 1, -1)
  )

  res <- RFGeneRank::factor_dependence(
    fit, se,
    covariates = "sex",
    method = "proxy",
    ngenes = 3,
    gene_selection = "importance"
  )

  expect_true(is.data.frame(res))
  expect_true(nrow(res) == 3)
  expect_true(all(res$covariate == "sex"))
})

test_that("factor_dependence errors when positive is invalid numeric", {
  expr <- matrix(rnorm(20), nrow = 5)
  rownames(expr) <- paste0("g", 1:5)
  colnames(expr) <- paste0("s", 1:4)

  cd <- S4Vectors::DataFrame(state = factor(c("A", "A", "B", "B")))
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(4), ncol = 1, dimnames = list(colnames(expr), "B"))
  fit <- methods::new("GeneRankFit", oof = list(prob = prob, y = cd$state))

  expect_error(
    RFGeneRank::factor_dependence(
      fit, se,
      covariates = "x",
      method = "proxy",
      gene_selection = "variance",
      ngenes = 2,
      positive = 3
    ),
    "`positive` must be NULL"
  )
})

test_that("factor_dependence errors when chosen positive label not in y levels", {
  expr <- matrix(rnorm(20), nrow = 5)
  rownames(expr) <- paste0("g", 1:5)
  colnames(expr) <- paste0("s", 1:4)

  cd <- S4Vectors::DataFrame(state = factor(c("A", "A", "B", "B")))
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(4), ncol = 1, dimnames = list(colnames(expr), "B"))
  fit <- methods::new("GeneRankFit", oof = list(prob = prob, y = cd$state))

  expect_error(
    RFGeneRank::factor_dependence(
      fit, se,
      covariates = "x",
      method = "proxy",
      gene_selection = "variance",
      ngenes = 2,
      positive = "Z"
    ),
    "Chosen positive class"
  )
})

test_that("factor_dependence drops state from covariates", {
  set.seed(1)

  expr <- matrix(rnorm(20 * 8), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:8)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 4)),
    sex = factor(rep(c("M", "F"), 4))
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(8), ncol = 1, dimnames = list(colnames(expr), "B"))
  fit <- methods::new("GeneRankFit", oof = list(prob = prob, y = cd$state))
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = seq_len(nrow(expr))
  )

  res <- RFGeneRank::factor_dependence(
    fit, se,
    covariates = c("state", "sex"),
    method = "proxy",
    ngenes = 5,
    gene_selection = "importance"
  )

  expect_true(is.data.frame(res))
  expect_false("state" %in% unique(res$covariate))
  expect_true("sex" %in% unique(res$covariate))
})

test_that("factor_dependence errors on missing expr assay", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("methods")

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(X = matrix(rnorm(20), nrow = 5)),
    colData = S4Vectors::DataFrame(state = factor(rep(c("A","B"), each = 2)))
  )

  prob <- matrix(runif(4), ncol = 1)
  rownames(prob) <- colnames(SummarizedExperiment::assay(se))
  colnames(prob) <- "B"

  fit <- methods::new("GeneRankFit", oof = list(prob = prob, y = factor(rep(c("A","B"), each = 2))))

  expect_error(
    RFGeneRank::factor_dependence(fit, se, covariates = "x", method = "proxy", ngenes = 2),
    "Assay 'expr' not found"
  )
})


test_that("factor_dependence errors on sample ID mismatch (assay expr)", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("methods")

  expr <- matrix(rnorm(5 * 4), nrow = 5)
  rownames(expr) <- paste0("g", 1:5)
  colnames(expr) <- paste0("s", 1:4)

  cd <- S4Vectors::DataFrame(state = factor(rep(c("A","B"), each = 2)))
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd)

  # oof contains an ID not in expr colnames
  prob <- matrix(runif(4), ncol = 1, dimnames = list(c("s1","s2","s3","BAD"), "B"))
  fit <- methods::new("GeneRankFit", oof = list(prob = prob, y = factor(rep(c("A","B"), each = 2))))

  expect_error(
    RFGeneRank::factor_dependence(fit, se, covariates = "x", method = "proxy", ngenes = 2),
    "Sample ID mismatch"
  )
})


test_that("factor_dependence errors when state missing or not binary", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("methods")

  expr <- matrix(rnorm(6 * 6), nrow = 6)
  rownames(expr) <- paste0("g", 1:6)
  colnames(expr) <- paste0("s", 1:6)

  cd1 <- S4Vectors::DataFrame(sex = factor(rep(c("M","F"), 3)))
  rownames(cd1) <- colnames(expr)
  se1 <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd1)

  prob <- matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B"))
  fit <- methods::new("GeneRankFit", oof = list(prob = prob, y = factor(rep(c("A","B"), 3))))

  expect_error(
    RFGeneRank::factor_dependence(fit, se1, covariates = "sex", method = "proxy", ngenes = 2),
    "colData\\(se\\)\\$state is required"
  )

  cd2 <- S4Vectors::DataFrame(state = factor(c("A","B","C","A","B","C")))
  rownames(cd2) <- colnames(expr)
  se2 <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd2)

  fit2 <- methods::new("GeneRankFit", oof = list(prob = prob, y = cd2$state))

  expect_error(
    RFGeneRank::factor_dependence(fit2, se2, covariates = "state", method = "proxy", ngenes = 2),
    "Binary classification expected"
  )
})


test_that("factor_dependence positive resolution branches error correctly", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("methods")

  expr <- matrix(rnorm(6 * 6), nrow = 6)
  rownames(expr) <- paste0("g", 1:6)
  colnames(expr) <- paste0("s", 1:6)

  cd <- S4Vectors::DataFrame(state = factor(rep(c("A","B"), each = 3)))
  rownames(cd) <- colnames(expr)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd)

  prob <- matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B"))
  fit <- methods::new("GeneRankFit", oof = list(prob = prob, y = cd$state))

  # positive numeric must be 1/2 only
  expect_error(
    RFGeneRank::factor_dependence(fit, se, covariates = "x", method = "proxy",
                                 positive = 3, ngenes = 2, gene_selection = "variance"),
    "`positive` must be NULL"
  )

  # positive label not in levels(y)
  expect_error(
    RFGeneRank::factor_dependence(fit, se, covariates = "x", method = "proxy",
                                 positive = "Z", ngenes = 2, gene_selection = "variance"),
    "Chosen positive class.*not in levels"
  )
})


test_that("factor_dependence importance/variance gene-selection branches", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("methods")

  expr <- matrix(rnorm(6 * 6), nrow = 6)
  rownames(expr) <- paste0("g", 1:6)
  colnames(expr) <- paste0("s", 1:6)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A","B"), each = 3)),
    sex   = factor(rep(c("M","F"), 3)),
    age   = c(10,20,30,40,50,60)
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd)

  prob <- matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B"))

  fit_bad_imp <- methods::new("GeneRankFit",
                              oof = list(prob = prob, y = cd$state),
                              imp = data.frame(x = 1))
  expect_error(
    RFGeneRank::factor_dependence(fit_bad_imp, se, covariates = "sex",
                                 method = "proxy", ngenes = 2, gene_selection = "importance"),
    "requires an importance table"
  )


  expr0 <- matrix(1, nrow = 6, ncol = 6)
  rownames(expr0) <- rownames(expr); colnames(expr0) <- colnames(expr)
  se0 <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr0), colData = cd)

  fit_any <- methods::new("GeneRankFit", oof = list(prob = prob, y = cd$state))

  expect_error(
    RFGeneRank::factor_dependence(fit_any, se0, covariates = "sex",
                                 method = "proxy", ngenes = 3, gene_selection = "variance"),
    "No variable genes found"
  )
})


test_that("factor_dependence covariate numeric vs factor and covariate filtering", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("methods")

  set.seed(1)

  expr <- matrix(rnorm(20 * 10), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:10)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A","B"), each = 5)),
    sex   = factor(rep(c("M","F"), length.out = 10)),
    age   = seq_len(10)
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = expr), colData = cd)

  prob <- matrix(runif(10), ncol = 1, dimnames = list(colnames(expr), "B"))

  imp <- data.frame(
    gene = rownames(expr),
    importance = seq_len(nrow(expr)),
    direction = rep(1, nrow(expr))
  )
  fit <- methods::new("GeneRankFit", oof = list(prob = prob, y = cd$state), imp = imp)

r0 <- RFGeneRank::factor_dependence(
  fit, se,
  covariates = c("state", "sex"),
  method = "proxy",
  ngenes = 5,
  gene_selection = "importance"
)

expect_true(is.data.frame(r0))
expect_false("state" %in% unique(r0$covariate))
expect_true("sex" %in% unique(r0$covariate))

  # 2) covariate factor (sex) + numeric (age) => 
  r1 <- RFGeneRank::factor_dependence(fit, se, covariates = c("sex","age"),
                                     method = "proxy", ngenes = 5, gene_selection = "importance")
  expect_true(is.data.frame(r1))
  expect_true(any(r1$test == "LM (v ~ y + covariate)"))
  expect_true(any(r1$test == "partial Spearman"))
})

toy_expr_samples_by_genes <- function(n_samples = 20, n_genes = 60, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n_samples * n_genes), nrow = n_samples, ncol = n_genes)
  colnames(X) <- paste0("g", seq_len(n_genes))
  rownames(X) <- paste0("s", seq_len(n_samples))
  X
}

toy_se_expr_genes_by_samples <- function(X_samp_by_gene, y, batch = NULL, age = NULL) {
  stopifnot(nrow(X_samp_by_gene) == length(y))
  if (is.null(batch)) batch <- factor(rep(c("b1","b2"), length.out = length(y)))
  if (is.null(age))   age   <- seq_len(length(y))

  expr <- t(X_samp_by_gene)  # genes x samples
  SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = y,
      batch = batch,
      age   = age,
      row.names = colnames(expr)
    )
  )
}

test_that("sign_importance(mean) returns signed_importance using fit+X+y", {
  skip_if_not_installed("methods")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  X <- toy_expr_samples_by_genes(n_samples = 24, n_genes = 80, seed = 1)
  y <- factor(rep(c("Control","Case"), each = 12))

  fit <- methods::new("GeneRankFit")

  # imp(fit) must be a data.frame with gene, importance (matches colnames(X))
  imp_df <- data.frame(
    gene = colnames(X),
    importance = seq_along(colnames(X)),   # deterministic, descending after sort
    stringsAsFactors = FALSE
  )
  imp(fit) <- imp_df

  out <- sign_importance(fit = fit, X = X, y = y, method = "mean", seed = 1)

  expect_true(is.data.frame(out))
  expect_true(all(c("gene","importance","direction","signed_importance") %in% names(out)))
  expect_true(nrow(out) > 0)
  expect_true(all(out$gene %in% colnames(X)))
})

test_that("sign_importance(de) works with de_table (gene, log2FC)", {
  skip_if_not_installed("methods")

  X <- toy_expr_samples_by_genes(n_samples = 20, n_genes = 40, seed = 2)
  y <- factor(rep(c("Control","Case"), each = 10))

  fit <- methods::new("GeneRankFit")
  imp(fit) <- data.frame(gene = colnames(X), importance = runif(ncol(X)))

  de_table <- data.frame(
    gene  = colnames(X)[1:20],
    log2FC = rnorm(20),
    stringsAsFactors = FALSE
  )

  out <- sign_importance(fit = fit, X = X, y = y, method = "de", de_table = de_table)

  expect_true(is.data.frame(out))
  expect_true(all(c("gene","log2FC","direction") %in% names(out)))
})

test_that("factor_dependence(method='proxy') returns results data.frame", {
  skip_if_not_installed("methods")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  # Build minimal fit with required oof$prob rownames == sample IDs
  n <- 30
  y <- factor(rep(c("Control","Case"), each = n/2))
  set.seed(3)
  p_case <- c(runif(n/2, 0.05, 0.35), runif(n/2, 0.65, 0.95))
  prob <- cbind(Control = 1 - p_case, Case = p_case)
  rownames(prob) <- paste0("s", seq_len(n))
  names(y) <- rownames(prob)

  fit <- methods::new("GeneRankFit")
  methods::slot(fit, "oof") <- list(prob = prob, y = y)

  # Provide importance table (required by proxy path)
  X <- toy_expr_samples_by_genes(n_samples = n, n_genes = 120, seed = 4)
  rownames(X) <- rownames(prob)  # align sample IDs

  imp(fit) <- data.frame(
    gene = colnames(X),
    importance = runif(ncol(X)),
    stringsAsFactors = FALSE
  )

  se <- toy_se_expr_genes_by_samples(
    X_samp_by_gene = X,
    y = y,
    batch = factor(rep(c("b1","b2","b3"), length.out = n)),
    age = sample(30:70, n, replace = TRUE)
  )

  res <- factor_dependence(
    fit = fit,
    se = se,
    covariates = c("batch", "age"),
    method = "proxy",
    ngenes = 50L,
    nsim = 16L,
    bg_per_class = 8L,
    seed = 1L
  )

  expect_true(is.data.frame(res))
  expect_true(all(c("covariate","gene","pval","fdr","dependent") %in% names(res)))
  expect_true(nrow(res) > 0)
  expect_true(attr(res, "method") %in% c("proxy","shap"))
})

test_that("factor_dependence keeps only requested existing covariates", {
  set.seed(10)

  expr <- matrix(rnorm(12 * 8), nrow = 12)
  rownames(expr) <- paste0("g", 1:12)
  colnames(expr) <- paste0("s", 1:8)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 4)),
    sex   = factor(rep(c("M", "F"), 4)),
    age   = c(21, 25, 29, 33, 37, 41, 45, 49)
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(8), ncol = 1, dimnames = list(colnames(expr), "B"))

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = cd$state)
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = seq_len(nrow(expr))
  )

  res <- RFGeneRank::factor_dependence(
    fit, se,
    covariates = c("sex", "age", "not_here"),
    method = "proxy",
    ngenes = 4,
    gene_selection = "importance"
  )

  expect_true(is.data.frame(res))
  expect_true(all(unique(res$covariate) %in% c("sex", "age")))
})

test_that("factor_dependence returns result with expected attributes for variance selection", {
  set.seed(11)

  expr <- matrix(rnorm(20 * 8), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:8)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 4)),
    age   = c(20, 24, 28, 32, 36, 40, 44, 48)
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(8), ncol = 1, dimnames = list(colnames(expr), "B"))

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = cd$state)
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = seq_len(nrow(expr))
  )

  res <- RFGeneRank::factor_dependence(
    fit, se,
    covariates = "age",
    method = "proxy",
    ngenes = 5,
    gene_selection = "variance"
  )

  expect_true(is.data.frame(res))
  expect_identical(attr(res, "method"), "proxy")
  expect_identical(attr(res, "assay"), "expr")
  expect_true(length(attr(res, "genes_used")) >= 1)
})

test_that("factor_dependence proxy path works with numeric positive index 2", {
  set.seed(12)

  expr <- matrix(rnorm(10 * 6), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:6)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 3)),
    age   = c(22, 26, 30, 34, 38, 42)
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B"))

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = cd$state)
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = rev(seq_len(nrow(expr)))
  )

  res <- RFGeneRank::factor_dependence(
    fit, se,
    covariates = "age",
    method = "proxy",
    ngenes = 3,
    gene_selection = "importance",
    positive = 2
  )

  expect_true(is.data.frame(res))
  expect_true(all(res$covariate == "age"))
})

test_that("factor_dependence returns non-tabular empty-like object when no usable covariates remain", {
  set.seed(13)

  expr <- matrix(rnorm(8 * 6), nrow = 8)
  rownames(expr) <- paste0("g", 1:8)
  colnames(expr) <- paste0("s", 1:6)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 3))
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(6), ncol = 1, dimnames = list(colnames(expr), "B"))

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = cd$state)
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = seq_len(nrow(expr))
  )

  res <- RFGeneRank::factor_dependence(
    fit, se,
    covariates = "state",
    method = "proxy",
    ngenes = 3,
    gene_selection = "importance"
  )

  expect_false(is.data.frame(res))
  expect_true(is.list(res) || is.null(res))
})

test_that("factor_dependence errors when top-importance genes do not overlap expr rownames", {
  expr <- matrix(rnorm(6 * 4), nrow = 6)
  rownames(expr) <- paste0("g", 1:6)
  colnames(expr) <- paste0("s", 1:4)

  cd <- S4Vectors::DataFrame(
    state = factor(c("A", "A", "B", "B")),
    sex   = factor(c("M", "F", "M", "F"))
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(4), ncol = 1, dimnames = list(colnames(expr), "B"))

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = cd$state)
  fit@imp <- data.frame(
    gene = c("x1", "x2", "x3"),
    importance = c(3, 2, 1)
  )

  expect_error(
    RFGeneRank::factor_dependence(
      fit, se,
      covariates = "sex",
      method = "proxy",
      ngenes = 2,
      gene_selection = "importance"
    ),
    "No overlap between top-importance genes and 'expr' rownames"
  )
})

test_that("factor_dependence errors when variance selection finds no variable genes", {
  expr <- matrix(1, nrow = 5, ncol = 4)
  rownames(expr) <- paste0("g", 1:5)
  colnames(expr) <- paste0("s", 1:4)

  cd <- S4Vectors::DataFrame(
    state = factor(c("A", "A", "B", "B")),
    age   = c(20, 30, 40, 50)
  )
  rownames(cd) <- colnames(expr)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  prob <- matrix(runif(4), ncol = 1, dimnames = list(colnames(expr), "B"))

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = cd$state)
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = seq_len(nrow(expr))
  )

  expect_error(
    RFGeneRank::factor_dependence(
      fit, se,
      covariates = "age",
      method = "proxy",
      ngenes = 3,
      gene_selection = "variance"
    ),
    "No variable genes found in assay 'expr'"
  )
})