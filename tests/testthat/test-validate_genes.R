test_that("validate_genes errors on missing genes argument", {
  expr <- matrix(rnorm(10 * 8), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:8)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = factor(rep(c("A", "B"), each = 4))
    )
  )

  expect_error(
    RFGeneRank::validate_genes(se = se, genes = character(0), methods = "ranger", k = 2),
    "non-empty.*genes"
  )
})

test_that("validate_genes errors when genes are missing from assay", {
  expr <- matrix(rnorm(10 * 8), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:8)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = factor(rep(c("A", "B"), each = 4))
    )
  )

  expect_error(
    RFGeneRank::validate_genes(se = se, genes = c("g1", "BAD"), methods = "ranger", k = 2),
    "Missing genes in assay\\(se\\)"
  )
})

test_that("validate_genes errors when label_col is not binary", {
  expr <- matrix(rnorm(10 * 9), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:9)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = factor(rep(c("A", "B", "C"), each = 3))
    )
  )

  expect_error(
    RFGeneRank::validate_genes(
      se = se,
      genes = rownames(expr)[1:5],
      methods = "ranger",
      k = 3,
      label_col = "state"
    ),
    "`label_col` must be a binary factor"
  )
})

test_that("validate_genes errors when positive is not in levels", {
  expr <- matrix(rnorm(10 * 8), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:8)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = factor(rep(c("A", "B"), each = 4))
    )
  )

  expect_error(
    RFGeneRank::validate_genes(
      se = se,
      genes = rownames(expr)[1:5],
      methods = "ranger",
      k = 2,
      label_col = "state",
      positive = "Z"
    ),
    "`positive` not found in levels\\(y\\)"
  )
})

test_that("validate_genes runs with ranger and returns summary plus method entry", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("caret")
  skip_if_not_installed("pROC")

  set.seed(1)

  expr <- matrix(rnorm(20 * 12), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:12)

  state <- factor(rep(c("A", "B"), each = 6))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = state)
  )

  res <- RFGeneRank::validate_genes(
    se = se,
    genes = rownames(expr)[1:6],
    methods = "ranger",
    k = 3,
    seed = 1,
    label_col = "state",
    calibrate = "none",
    thr_metric = "youden"
  )

  expect_true(is.list(res))
  expect_true("summary" %in% names(res))
  expect_true("ranger" %in% names(res))

  expect_true(is.data.frame(res$summary))
  expect_true(is.list(res$ranger))
  expect_true(all(c("oof", "threshold", "conf_mat", "metrics") %in% names(res$ranger)))
  expect_true(is.data.frame(res$ranger$oof))
  expect_true(all(c("p_raw", "p_cal", "p_use", "y", "fold") %in% names(res$ranger$oof)))
})

test_that("validate_genes runs with positive relevel branch", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("caret")
  skip_if_not_installed("pROC")

  set.seed(2)

  expr <- matrix(rnorm(20 * 12), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:12)

  state <- factor(rep(c("A", "B"), each = 6))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = state)
  )

  res <- RFGeneRank::validate_genes(
    se = se,
    genes = rownames(expr)[1:6],
    methods = "ranger",
    k = 3,
    seed = 1,
    label_col = "state",
    positive = "A",
    calibrate = "none"
  )

  expect_true(is.list(res))
  expect_true("ranger" %in% names(res))
  expect_true(is.data.frame(res$ranger$oof))
})

test_that("validate_genes runs cost threshold branch", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("caret")
  skip_if_not_installed("pROC")

  set.seed(3)

  expr <- matrix(rnorm(20 * 12), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:12)

  state <- factor(rep(c("A", "B"), each = 6))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = state)
  )

  res <- RFGeneRank::validate_genes(
    se = se,
    genes = rownames(expr)[1:6],
    methods = "ranger",
    k = 3,
    seed = 1,
    label_col = "state",
    calibrate = "none",
    thr_metric = "cost",
    cost = c(fp = 2, fn = 1)
  )

  expect_true(is.list(res))
  expect_true(is.numeric(res$ranger$threshold))
  expect_equal(length(res$ranger$threshold), 1)
})

test_that("validate_genes errors on empty genes and missing genes", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  expr <- matrix(rnorm(10 * 8), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:8)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = factor(rep(c("A", "B"), each = 4)))
  )

 expect_error(
  RFGeneRank::validate_genes(se = se, genes = character(0), methods = "ranger", k = 2),
  "non-empty character vector genes"
)

  expect_error(
    RFGeneRank::validate_genes(se = se, genes = c("g1", "BAD"), methods = "ranger", k = 2),
    "Missing genes in assay\\(se\\)"
  )
})

test_that("validate_genes errors when label_col is not binary or positive is absent", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  expr <- matrix(rnorm(10 * 9), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:9)

  se3 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = factor(rep(c("A", "B", "C"), each = 3)))
  )

  expect_error(
    RFGeneRank::validate_genes(
      se = se3,
      genes = rownames(expr)[1:5],
      methods = "ranger",
      k = 3,
      label_col = "state"
    ),
    "`label_col` must be a binary factor"
  )

  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr[, 1:8, drop = FALSE]),
    colData = S4Vectors::DataFrame(state = factor(rep(c("A", "B"), each = 4)))
  )

  expect_error(
    RFGeneRank::validate_genes(
      se = se2,
      genes = rownames(expr)[1:5],
      methods = "ranger",
      k = 2,
      label_col = "state",
      positive = "Z"
    ),
    "`positive` not found in levels\\(y\\)"
  )
})

test_that("validate_genes errors on unsupported methods", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  expr <- matrix(rnorm(10 * 8), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:8)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = factor(rep(c("A", "B"), each = 4)))
  )

  expect_error(
    RFGeneRank::validate_genes(
      se = se,
      genes = rownames(expr)[1:5],
      methods = c("ranger", "badmethod"),
      k = 2,
      label_col = "state"
    ),
    "Unsupported methods"
  )
})

test_that("validate_genes runs ranger path and returns expected components", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  skip_if_not_installed("ranger")
  skip_if_not_installed("caret")
  skip_if_not_installed("pROC")

  set.seed(1)

  expr <- matrix(rnorm(20 * 12), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:12)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = factor(rep(c("A", "B"), each = 6)))
  )

  res <- RFGeneRank::validate_genes(
    se = se,
    genes = rownames(expr)[1:6],
    methods = "ranger",
    k = 3,
    seed = 1,
    label_col = "state",
    calibrate = "none",
    thr_metric = "youden"
  )

  expect_true(is.list(res))
  expect_true("summary" %in% names(res))
  expect_true("ranger" %in% names(res))
  expect_true(is.data.frame(res$summary))
  expect_true(is.list(res$ranger))
  expect_true(all(c("oof", "threshold", "conf_mat", "metrics") %in% names(res$ranger)))
  expect_true(all(c("p_raw", "p_cal", "p_use", "y", "fold") %in% names(res$ranger$oof)))
})

test_that("validate_genes runs positive relevel and cost-threshold branches", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  skip_if_not_installed("ranger")
  skip_if_not_installed("caret")
  skip_if_not_installed("pROC")

  set.seed(2)

  expr <- matrix(rnorm(20 * 12), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:12)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = factor(rep(c("A", "B"), each = 6)))
  )

  res <- RFGeneRank::validate_genes(
    se = se,
    genes = rownames(expr)[1:6],
    methods = "ranger",
    k = 3,
    seed = 1,
    label_col = "state",
    positive = "A",
    calibrate = "none",
    thr_metric = "cost",
    cost = c(fp = 2, fn = 1)
  )

  expect_true(is.list(res))
  expect_true(is.numeric(res$ranger$threshold))
  expect_equal(length(res$ranger$threshold), 1)
})

test_that("validate_genes f1 threshold branch currently errors through pROC::coords", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  skip_if_not_installed("ranger")
  skip_if_not_installed("caret")
  skip_if_not_installed("pROC")

  set.seed(3)

  expr <- matrix(rnorm(20 * 12), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:12)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = factor(rep(c("A", "B"), each = 6)))
  )

  expect_error(
    RFGeneRank::validate_genes(
      se = se,
      genes = rownames(expr)[1:6],
      methods = "ranger",
      k = 3,
      seed = 1,
      label_col = "state",
      calibrate = "none",
      thr_metric = "f1"
    ),
    "'arg' should be one of"
  )
})

test_that("validate_genes supports repeated methods input via unique()", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  skip_if_not_installed("ranger")
  skip_if_not_installed("caret")
  skip_if_not_installed("pROC")

  set.seed(4)

  expr <- matrix(rnorm(20 * 12), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:12)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = factor(rep(c("A", "B"), each = 6)))
  )

  res <- RFGeneRank::validate_genes(
    se = se,
    genes = rownames(expr)[1:6],
    methods = c("ranger", "ranger"),
    k = 3,
    seed = 1,
    label_col = "state",
    calibrate = "none"
  )

  expect_true(is.list(res))
  expect_equal(sum(names(res) == "ranger"), 1)
})

test_that("validate_genes runs with ranger-only on larger balanced toy SE", {
  skip_if_not_installed("SummarizedExperiment")

  se <- toy_se_cov(ngenes = 200, nsamp = 80, seed = 20, entrez = FALSE)

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  cd$batch <- factor(rep(c("B1","B2"), each = nrow(cd)/2))
  cd$state <- factor(rep(c("Control","Case"), length.out = nrow(cd)))
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cd)

  fit <- RFGeneRank::rank_genes(
    se = se,
    label_col = "state",
    batch_col = "batch",
    trees = 50,
    n_top = 80,
    cv = "kfold"
  )

  vg <- RFGeneRank::validate_genes
  f <- names(formals(vg))

  # gene set
  gset <- rownames(SummarizedExperiment::assay(se, "expr"))[1:40]

  args <- list()
  if ("fit" %in% f) args$fit <- fit
  if ("se" %in% f) args$se <- se
  if ("genes" %in% f) args$genes <- gset
  if ("gene_set" %in% f) args$gene_set <- gset
  if ("label_col" %in% f) args$label_col <- "state"
  if ("batch_col" %in% f) args$batch_col <- "batch"
  if ("k" %in% f) args$k <- 3
  if ("seed" %in% f) args$seed <- 1

  # ranger-only
  if ("methods" %in% f) args$methods <- "ranger"
  if ("models" %in% f) args$models <- "ranger"
  if ("method" %in% f) args$method <- "ranger"
  if ("model" %in% f) args$model <- "ranger"

  # threshold 
  if ("threshold" %in% f) args$threshold <- 0.5
  if ("thr" %in% f) args$thr <- 0.5
  if ("thr_metric" %in% f) args$thr_metric <- "youden"
  if ("cost" %in% f) args$cost <- c(fp = 1, fn = 1)

  out <- do.call(vg, args)

  expect_true(is.list(out) || is.data.frame(out))
})

test_that("validate_genes rejects unknown methods early", {
  se <- toy_se(ngenes = 30, nsamp = 12, seed = 2)
  genes <- rownames(SummarizedExperiment::assay(se, "expr"))[1:10]

  expect_error(
    validate_genes(
      se = se,
      genes = genes,
      methods = "UNKNOWN_METHOD",
      k = 3,
      seed = 1,
      label_col = "label"
    ),
    "Unsupported methods"
  )
})