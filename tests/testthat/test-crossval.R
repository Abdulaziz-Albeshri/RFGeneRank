test_that("rfgr_crossval runs on toy data and returns expected keys", {
  skip_if_not_installed("SummarizedExperiment")

  se <- toy_se(ngenes = 60, nsamp = 20, seed = 1)
  expr <- SummarizedExperiment::assay(se, "expr")
  md <- as.data.frame(SummarizedExperiment::colData(se))

  out <- RFGeneRank::rfgr_crossval(
    expr = expr,
    metadata = md,
    label_col = "label",
    batch_col = "batch",
    covariates = NULL,
    cv = "kfold",
    k = 3,
    combat_mode = "none",
    rf_trees = 50,
    seed = 1,
    verbose = FALSE
  )

  expect_true(is.list(out))
  expect_true(all(c("auc_by_fold","mean_auc","settings","folds_info") %in% names(out)))
})

test_that("rfgr_crossval runs with kfold and returns expected structure", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")

  set.seed(1)

  expr <- matrix(rnorm(20 * 8), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:8)

  metadata <- data.frame(
    state = factor(rep(c("A", "B"), each = 4)),
    batch = factor(rep(c("b1", "b2"), times = 4)),
    row.names = colnames(expr)
  )

  res <- RFGeneRank::rfgr_crossval(
    expr = expr,
    metadata = metadata,
    label_col = "state",
    batch_col = "batch",
    cv = "kfold",
    k = 2,
    combat_mode = "none",
    rf_trees = 10,
    verbose = FALSE
  )

  expect_true(is.list(res))
  expect_true(all(c("auc_by_fold", "mean_auc", "settings", "folds_info") %in% names(res)))
  expect_true(is.numeric(res$auc_by_fold))
  expect_equal(length(res$auc_by_fold), 2)
  expect_true(is.numeric(res$mean_auc))
  expect_equal(res$settings$cv, "kfold")
  expect_equal(res$settings$combat_mode, "none")
  expect_equal(length(res$folds_info), 2)
})

test_that("rfgr_crossval runs with lobo and returns one fold per batch", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")

  set.seed(2)

  expr <- matrix(rnorm(20 * 8), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:8)

  metadata <- data.frame(
    state = factor(c("A", "A", "B", "B", "A", "A", "B", "B")),
    batch = factor(c("b1", "b1", "b1", "b1", "b2", "b2", "b2", "b2")),
    row.names = colnames(expr)
  )

  res <- RFGeneRank::rfgr_crossval(
    expr = expr,
    metadata = metadata,
    label_col = "state",
    batch_col = "batch",
    cv = "lobo",
    combat_mode = "none",
    rf_trees = 10,
    verbose = FALSE
  )

  expect_true(is.numeric(res$auc_by_fold))
  expect_equal(length(res$auc_by_fold), 2)
  expect_true(all(grepl("^test=", names(res$auc_by_fold))))
  expect_equal(length(res$folds_info), 2)
})

test_that("rfgr_crossval runs with groupk and returns expected fold count", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")

  set.seed(3)

  expr <- matrix(rnorm(20 * 12), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:12)

  metadata <- data.frame(
    state = factor(rep(c("A", "B"), each = 6)),
    batch = factor(rep(c("b1", "b2", "b3"), each = 4)),
    row.names = colnames(expr)
  )

  res <- RFGeneRank::rfgr_crossval(
    expr = expr,
    metadata = metadata,
    label_col = "state",
    batch_col = "batch",
    cv = "groupk",
    k = 2,
    combat_mode = "none",
    rf_trees = 10,
    verbose = FALSE
  )

  expect_true(is.numeric(res$auc_by_fold))
  expect_equal(length(res$auc_by_fold), 2)
  expect_true(all(grepl("^fold=", names(res$auc_by_fold))))
  expect_equal(res$settings$cv, "groupk")
})

test_that("rfgr_crossval accepts kfold without batch column", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("pROC")

  set.seed(4)

  expr <- matrix(rnorm(20 * 8), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:8)

  metadata <- data.frame(
    state = factor(rep(c("A", "B"), each = 4)),
    row.names = colnames(expr)
  )

  res <- RFGeneRank::rfgr_crossval(
    expr = expr,
    metadata = metadata,
    label_col = "state",
    cv = "kfold",
    k = 2,
    combat_mode = "none",
    rf_trees = 10,
    verbose = FALSE
  )

  expect_true(is.list(res))
  expect_equal(res$settings$cv, "kfold")
})

test_that("rfgr_crossval errors when batch_col is missing for lobo/groupk", {
  expr <- matrix(rnorm(10 * 6), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:6)

  metadata <- data.frame(
    state = factor(rep(c("A", "B"), each = 3)),
    row.names = colnames(expr)
  )

  expect_error(
    RFGeneRank::rfgr_crossval(
      expr = expr,
      metadata = metadata,
      label_col = "state",
      batch_col = "batch",
      cv = "lobo",
      combat_mode = "none",
      rf_trees = 10,
      verbose = FALSE
    )
  )

  expect_error(
    RFGeneRank::rfgr_crossval(
      expr = expr,
      metadata = metadata,
      label_col = "state",
      batch_col = "batch",
      cv = "groupk",
      k = 2,
      combat_mode = "none",
      rf_trees = 10,
      verbose = FALSE
    )
  )
})

test_that("rfgr_crossval errors when label_col is missing", {
  expr <- matrix(rnorm(10 * 6), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:6)

  metadata <- data.frame(
    batch = factor(rep(c("b1", "b2"), each = 3)),
    row.names = colnames(expr)
  )

  expect_error(
    RFGeneRank::rfgr_crossval(
      expr = expr,
      metadata = metadata,
      label_col = "state",
      batch_col = "batch",
      cv = "kfold",
      k = 2,
      combat_mode = "none",
      rf_trees = 10,
      verbose = FALSE
    )
  )
})

test_that("rfgr_crossval errors when expr and metadata dimensions mismatch", {
  expr <- matrix(rnorm(10 * 6), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:6)

  metadata <- data.frame(
    state = factor(rep(c("A", "B"), c(2, 3))),
    batch = factor(rep(c("b1", "b2"), length.out = 5)),
    row.names = paste0("s", 1:5)
  )

  expect_error(
    RFGeneRank::rfgr_crossval(
      expr = expr,
      metadata = metadata,
      label_col = "state",
      batch_col = "batch",
      cv = "kfold",
      k = 2,
      combat_mode = "none",
      rf_trees = 10,
      verbose = FALSE
    )
  )
})

.rfgr_call_crossval <- function(se, label_col = "state", batch_col = "batch", seed = 1, ...) {
  stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))

  fn <- RFGeneRank::rfgr_crossval
  f  <- names(formals(fn))

  expr <- SummarizedExperiment::assay(se, "expr")            # genes x samples
  meta <- as.data.frame(SummarizedExperiment::colData(se))   # samples x covariates
  y    <- SummarizedExperiment::colData(se)[[label_col]]

  # ensure sample IDs consistent (important for strict checks)
  if (!is.null(colnames(expr)) && is.null(rownames(meta))) {
    rownames(meta) <- colnames(expr)
  }
  if (!is.null(colnames(expr)) && !is.null(rownames(meta))) {
    meta <- meta[colnames(expr), , drop = FALSE]
  }

  args <- list(...)

  # --- Provide data in the form supported by rfgr_crossval() ---
  # Primary interface in your version appears to be: expr + metadata
  if ("expr" %in% f) args$expr <- expr
  if ("expr_mat" %in% f) args$expr_mat <- expr
  if ("metadata" %in% f) args$metadata <- meta
  if ("meta" %in% f) args$meta <- meta

  # Some versions also accept X/y
  if ("X" %in% f) args$X <- t(expr)  # samples x genes
  if ("y" %in% f) args$y <- as.factor(y)

  # Column-name parameters (only if supported)
  if ("label_col" %in% f) args$label_col <- label_col
  if ("state_col" %in% f) args$state_col <- label_col
  if ("outcome" %in% f) args$outcome <- label_col

  if ("batch_col" %in% f) args$batch_col <- batch_col
  if ("batch" %in% f) args$batch <- batch_col

  if ("seed" %in% f) args$seed <- seed

  # CV controls (only if supported)
  if ("cv" %in% f) args$cv <- "kfold"
  if ("k" %in% f) args$k <- 5
  if ("folds" %in% f) args$folds <- 5

  # Trees (only if supported)
  if ("trees" %in% f) args$trees <- 50
  if ("num.trees" %in% f) args[["num.trees"]] <- 50
  if ("num_trees" %in% f) args[["num_trees"]] <- 50
  if ("ntree" %in% f) args$ntree <- 50

  # Keep only args that exist in signature (prevents "unused arguments")
  args <- args[names(args) %in% f]

  do.call(fn, args)
}