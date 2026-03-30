test_that("prepare_data returns SummarizedExperiment", {
  se <- toy_se(ngenes = 60, nsamp = 12, seed = 1)

  se2 <- prepare_data(
    mats = list(SummarizedExperiment::assay(se, "expr")),
    metas = list(as.data.frame(SummarizedExperiment::colData(se))),
    label_col = "label",
    batch_col = "batch",
    n_var = 30,
    log1p = FALSE,
    batch_method = "none",
    filter_in_cv = TRUE
  )

  expect_true(inherits(se2, "SummarizedExperiment"))
  expect_true("expr" %in% SummarizedExperiment::assayNames(se2))
})

test_that("prepare_data applies log1p and n_var filter when filter_in_cv is FALSE", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  skip_if_not_installed("matrixStats")

  set.seed(1)
  expr <- abs(matrix(rnorm(100 * 12), nrow = 100))
  rownames(expr) <- paste0("g", 1:100)
  colnames(expr) <- paste0("s", 1:12)

  meta <- data.frame(
    label = factor(rep(c("A", "B"), each = 6)),
    batch = factor(rep(c("b1", "b2"), each = 6)),
    row.names = colnames(expr)
  )

  se <- suppressWarnings(
    RFGeneRank::prepare_data(
      mats = list(expr),
      metas = list(meta),
      label_col = "label",
      batch_col = "batch",
      n_var = 25,
      log1p = TRUE,
      batch_method = "none",
      filter_in_cv = FALSE
    )
  )

  x <- SummarizedExperiment::assay(se, "expr")
  expect_equal(nrow(x), 25)
  expect_true(all(is.finite(x)))
  expect_true(all(x >= 0))
})

test_that("prepare_data skips global variance filter when filter_in_cv is TRUE", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  set.seed(2)
  expr <- abs(matrix(rnorm(80 * 12), nrow = 80))
  rownames(expr) <- paste0("g", 1:80)
  colnames(expr) <- paste0("s", 1:12)

  meta <- data.frame(
    label = factor(rep(c("A", "B"), each = 6)),
    batch = factor(rep(c("b1", "b2"), each = 6)),
    row.names = colnames(expr)
  )

  se <- suppressWarnings(
    RFGeneRank::prepare_data(
      mats = list(expr),
      metas = list(meta),
      label_col = "label",
      batch_col = "batch",
      n_var = 20,
      log1p = FALSE,
      batch_method = "none",
      filter_in_cv = TRUE
    )
  )

  x <- SummarizedExperiment::assay(se, "expr")
  expect_equal(nrow(x), 80)
})

test_that("prepare_data stores metadata flags for no correction", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  expr <- matrix(rnorm(20 * 12), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:12)

  meta <- data.frame(
    label = factor(rep(c("A", "B"), each = 6)),
    batch = factor(rep(c("b1", "b2"), each = 6)),
    row.names = colnames(expr)
  )

  se <- suppressWarnings(
    RFGeneRank::prepare_data(
      mats = list(expr),
      metas = list(meta),
      label_col = "label",
      batch_col = "batch",
      n_var = 0,
      log1p = FALSE,
      batch_method = "none",
      batch_correction_scope = "global"
    )
  )

  md <- S4Vectors::metadata(se)
  expect_false(isTRUE(md$`..batch_corrected`))
  expect_identical(md$`..batch_correction_scope`, "none")
  expect_identical(md$`..batch_col`, "batch")
})

test_that("prepare_data stores fold scope when requested", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  expr <- matrix(rnorm(20 * 12), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:12)

  meta <- data.frame(
    label = factor(rep(c("A", "B"), each = 6)),
    batch = factor(rep(c("b1", "b2"), each = 6)),
    row.names = colnames(expr)
  )

  se <- suppressWarnings(
    RFGeneRank::prepare_data(
      mats = list(expr),
      metas = list(meta),
      label_col = "label",
      batch_col = "batch",
      n_var = 0,
      log1p = FALSE,
      batch_method = "none",
      batch_correction_scope = "fold",
      batch_covariates = c("label")
    )
  )

  md <- S4Vectors::metadata(se)
  expect_false(isTRUE(md$`..batch_corrected`))
  expect_identical(md$`..batch_correction_scope`, "fold")
  expect_identical(md$`..batch_covariates`, "label")
})

test_that("prepare_data warns on confounding diagnostics", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  expr <- matrix(rnorm(30 * 8), nrow = 30)
  rownames(expr) <- paste0("g", 1:30)
  colnames(expr) <- paste0("s", 1:8)

  meta <- data.frame(
    label = factor(c(rep("A", 5), rep("B", 3))),
    batch = factor(c(rep("b1", 5), rep("b2", 3))),
    row.names = colnames(expr)
  )

  expect_warning(
    RFGeneRank::prepare_data(
      mats = list(expr),
      metas = list(meta),
      label_col = "label",
      batch_col = "batch",
      n_var = 0,
      log1p = FALSE,
      batch_method = "none"
    ),
    "Strong label-batch association|Some batches have <5 samples|Empty cells|Chi-squared approximation may be incorrect"
  )
})

test_that("prepare_data errors on invalid batch_method and invalid scope", {
  expr <- matrix(rnorm(20 * 8), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:8)

  meta <- data.frame(
    label = factor(rep(c("A", "B"), each = 4)),
    batch = factor(rep(c("b1", "b2"), each = 4)),
    row.names = colnames(expr)
  )

  expect_error(
    RFGeneRank::prepare_data(
      mats = list(expr),
      metas = list(meta),
      label_col = "label",
      batch_col = "batch",
      batch_method = "bad"
    )
  )

  expect_error(
    RFGeneRank::prepare_data(
      mats = list(expr),
      metas = list(meta),
      label_col = "label",
      batch_col = "batch",
      batch_correction_scope = "bad"
    )
  )
})

test_that("prepare_data sets metadata flags for no correction", {
  expr <- matrix(rnorm(20 * 12), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:12)

  meta <- data.frame(
    label = factor(rep(c("A", "B"), each = 6)),
    batch = factor(rep(c("b1", "b2"), each = 6)),
    row.names = colnames(expr)
  )

  se <- RFGeneRank::prepare_data(
    mats = list(expr),
    metas = list(meta),
    label_col = "label",
    batch_col = "batch",
    n_var = 0,
    log1p = FALSE,
    batch_method = "none"
  )

  md <- S4Vectors::metadata(se)
  expect_false(isTRUE(md$`..batch_corrected`))
  expect_true(md$`..batch_correction_scope` %in% c("none", "fold"))
})

test_that("prepare_data warns on confounding diagnostics", {
  expr <- matrix(rnorm(20 * 8), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:8)

  meta <- data.frame(
    label = factor(c(rep("A", 5), rep("B", 3))),
    batch = factor(c(rep("b1", 5), rep("b2", 3))),
    row.names = colnames(expr)
  )

  expect_warning(
    RFGeneRank::prepare_data(
      mats = list(expr),
      metas = list(meta),
      label_col = "label",
      batch_col = "batch",
      n_var = 0,
      log1p = FALSE,
      batch_method = "none"
    ),
    "Strong label-batch association|Some batches have <5 samples|Empty cells"
  )
})

test_that("prepare_data applies log1p and n_var filtering branches", {
  expr <- matrix(abs(rnorm(100 * 12)), nrow = 100)
  rownames(expr) <- paste0("g", 1:100)
  colnames(expr) <- paste0("s", 1:12)

  meta <- data.frame(
    label = factor(rep(c("A", "B"), each = 6)),
    batch = factor(rep(c("b1", "b2"), each = 6)),
    row.names = colnames(expr)
  )

  se <- RFGeneRank::prepare_data(
    mats = list(expr),
    metas = list(meta),
    label_col = "label",
    batch_col = "batch",
    n_var = 25,
    log1p = TRUE,
    batch_method = "none"
  )

  x <- SummarizedExperiment::assay(se, "expr")
  expect_equal(nrow(x), 25)
  expect_true(all(is.finite(x)))
})

test_that("prepare_data works with correctly aligned metadata", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  skip_if_not_installed("matrixStats")

  expr <- matrix(rnorm(20 * 10), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:10)

  meta <- data.frame(
    label = factor(rep(c("A", "B"), each = 5)),
    batch = factor(rep(c("b1", "b2"), times = 5)),
    row.names = colnames(expr)
  )

  se <- RFGeneRank::prepare_data(
    mats = list(expr),
    metas = list(meta),
    label_col = "label",
    batch_col = "batch",
    n_var = 0,
    log1p = FALSE,
    batch_method = "none"
  )

  expect_true(inherits(se, "SummarizedExperiment"))
  expect_equal(colnames(SummarizedExperiment::assay(se, "expr")), rownames(as.data.frame(SummarizedExperiment::colData(se))))
})

test_that("prepare_data emits confounding warnings when conditions are met", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  skip_if_not_installed("matrixStats")

  expr <- matrix(rnorm(30 * 8), nrow = 30)
  rownames(expr) <- paste0("g", 1:30)
  colnames(expr) <- paste0("s", 1:8)

  meta <- data.frame(
    label = factor(c(rep("A", 5), rep("B", 3))),
    batch = factor(c(rep("b1", 5), rep("b2", 3))),
    row.names = colnames(expr)
  )

  expect_warning(
    se <- RFGeneRank::prepare_data(
      mats = list(expr),
      metas = list(meta),
      label_col = "label",
      batch_col = "batch",
      n_var = 0,
      log1p = FALSE,
      batch_method = "none"
    ),
    "Strong label-batch association|Some batches have <5 samples|Empty cells"
  )

  md <- S4Vectors::metadata(se)
  expect_false(isTRUE(md$`..batch_corrected`))
})