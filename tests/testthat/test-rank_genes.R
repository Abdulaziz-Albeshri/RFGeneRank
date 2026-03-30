test_that("rank_genes returns GeneRankFit and key outputs are accessible", {
  se <- toy_se(ngenes = 40, nsamp = 12, seed = 1)

  fit <- rank_genes(
    se = se,
    label_col = "label",
    n_top = 15,
    trees = 50,
    seed = 1
  )

  expect_s4_class(fit, "GeneRankFit")

  # accessors should not error
  imptab <- imp(fit)
  expect_true(is.null(imptab) || is.data.frame(imptab))

  o <- oof(fit)
  expect_true(is.null(o) || is.list(o))

  p <- params(fit)
  expect_true(is.list(p))

  # top_genes: allow character, data.frame, or list
  tg <- top_genes(fit, n = 10, map = FALSE)

  expect_true(
    is.character(tg) || is.data.frame(tg) || is.list(tg),
    info = paste("top_genes returned class:", paste(class(tg), collapse = ", "))
  )

  if (is.character(tg)) {
    expect_true(length(tg) <= 10)
    expect_true(all(nzchar(tg)))
  } else if (is.data.frame(tg)) {
    expect_true(nrow(tg) <= 10)
    expect_true(any(c("gene", "Gene", "id") %in% colnames(tg)))
  } else if (is.list(tg)) {
    # try to detect a gene vector inside the list
    g <- NULL
    for (nm in c("genes", "gene", "top_genes", "ids", "features")) {
      if (!is.null(tg[[nm]])) {
        g <- tg[[nm]]
        break
      }
    }
    if (is.null(g) && length(tg) == 1 && is.character(tg[[1]])) {
      g <- tg[[1]]
    }

    expect_true(!is.null(g), info = "top_genes returned a list but no recognizable gene vector found.")
    expect_true(is.character(g))
    expect_true(length(g) > 0)
  }
})

test_that("rank_genes errors when label_col is missing", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ranger")

  expr <- matrix(rnorm(20 * 8), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:8)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(group = factor(rep(c("A", "B"), each = 4)))
  )

  expect_error(
    RFGeneRank::rank_genes(se, label_col = "state", trees = 10, k = 2, auto_confounds = FALSE),
    "label_col not found in colData\\(se\\)"
  )
})

test_that("rank_genes errors when there are fewer than 2 classes", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ranger")

  expr <- matrix(rnorm(20 * 6), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:6)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = factor(rep("A", 6)))
  )

  expect_error(
    RFGeneRank::rank_genes(se, label_col = "state", trees = 10, k = 2, auto_confounds = FALSE),
    "Need at least 2 classes"
  )
})

test_that("rank_genes blocks double batch correction when prepare_data metadata says global correction already happened", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ranger")

  expr <- matrix(rnorm(20 * 8), nrow = 20)
  rownames(expr) <- paste0("g", 1:20)
  colnames(expr) <- paste0("s", 1:8)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = factor(rep(c("A", "B"), each = 4)),
      batch = factor(rep(c("b1", "b2"), 4))
    )
  )
  S4Vectors::metadata(se)$`..batch_corrected` <- TRUE
  S4Vectors::metadata(se)$`..batch_correction_scope` <- "global"

  expect_error(
    RFGeneRank::rank_genes(
      se,
      label_col = "state",
      batch_col = "batch",
      fold_batch_correction = TRUE,
      trees = 10,
      k = 2,
      auto_confounds = FALSE
    ),
    "Batch correction appears to have been applied in prepare_data"
  )
})

test_that("rank_genes runs and returns GeneRankFit with expected slots", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ranger")

  set.seed(1)

  expr <- matrix(rnorm(40 * 8), nrow = 40)
  rownames(expr) <- paste0("g", 1:40)
  colnames(expr) <- paste0("s", 1:8)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = factor(rep(c("A", "B"), each = 4))
    )
  )

  fit <- RFGeneRank::rank_genes(
    se,
    label_col = "state",
    trees = 10,
    k = 2,
    n_top = 10,
    auto_confounds = FALSE,
    cv = "kfold"
  )

  expect_s4_class(fit, "GeneRankFit")
  expect_true(is.list(fit@params))
  expect_true(is.list(fit@oof))
  expect_true(is.data.frame(fit@imp))
  expect_true(is.character(fit@features))
  expect_true(length(fit@features) > 0)
  expect_true(all(c("prob", "pred", "y") %in% names(fit@oof)))
  expect_equal(nrow(fit@oof$prob), ncol(expr))
  expect_equal(colnames(fit@oof$prob), levels(factor(SummarizedExperiment::colData(se)$state)))
  expect_true(all(c("gene", "importance", "SelectedInFolds") %in% names(fit@imp)))
})

test_that("rank_genes honors transform=log1p and standardize branches", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ranger")
  skip_if_not_installed("matrixStats")

  set.seed(2)

  expr <- abs(matrix(rnorm(40 * 8), nrow = 40))
  rownames(expr) <- paste0("g", 1:40)
  colnames(expr) <- paste0("s", 1:8)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = factor(rep(c("A", "B"), each = 4))
    )
  )

  fit <- RFGeneRank::rank_genes(
    se,
    label_col = "state",
    trees = 10,
    k = 2,
    n_top = 12,
    transform = "log1p",
    standardize = TRUE,
    auto_confounds = FALSE,
    cv = "kfold"
  )

  expect_s4_class(fit, "GeneRankFit")
  expect_identical(fit@params$transform, "log1p")
  expect_true(isTRUE(fit@params$standardize))
})

test_that("rank_genes honors filter_low_expr branch", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ranger")

  set.seed(3)

  expr <- matrix(0, nrow = 30, ncol = 8)
  expr[1:10, ] <- matrix(rpois(10 * 8, lambda = 5), nrow = 10)
  rownames(expr) <- paste0("g", 1:30)
  colnames(expr) <- paste0("s", 1:8)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = factor(rep(c("A", "B"), each = 4))
    )
  )

  fit <- RFGeneRank::rank_genes(
    se,
    label_col = "state",
    trees = 10,
    k = 2,
    filter_low_expr = TRUE,
    min_prop = 0.5,
    n_top = 5,
    auto_confounds = FALSE,
    cv = "kfold"
  )

  expect_s4_class(fit, "GeneRankFit")
  expect_true(isTRUE(fit@params$filter_low_expr))
  expect_equal(fit@params$min_prop, 0.5)
})

test_that("rank_genes computes class_weights automatically when imbalance is >= 1.5x", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ranger")

  set.seed(4)

  expr <- matrix(rnorm(40 * 10), nrow = 40)
  rownames(expr) <- paste0("g", 1:40)
  colnames(expr) <- paste0("s", 1:10)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = factor(c(rep("A", 7), rep("B", 3)))
    )
  )

  fit <- RFGeneRank::rank_genes(
    se,
    label_col = "state",
    trees = 10,
    k = 2,
    n_top = 10,
    auto_confounds = FALSE,
    cv = "kfold"
  )

  expect_s4_class(fit, "GeneRankFit")
  expect_false(is.null(fit@params$class_weights))
  expect_true(is.numeric(fit@params$class_weights))
})

test_that("rank_genes recovers batch defaults from colData metadata columns", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ranger")

  set.seed(5)

  expr <- matrix(rnorm(40 * 8), nrow = 40)
  rownames(expr) <- paste0("g", 1:40)
  colnames(expr) <- paste0("s", 1:8)

  cd <- S4Vectors::DataFrame(
    state = factor(rep(c("A", "B"), each = 4)),
    batch = factor(rep(c("b1", "b2"), 4))
  )
  cd$..batch_col <- rep("batch", 8)
  cd$..batch_covariates <- I(replicate(8, c("state"), simplify = FALSE))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )

  fit <- RFGeneRank::rank_genes(
    se,
    label_col = "state",
    trees = 10,
    k = 2,
    n_top = 10,
    auto_confounds = FALSE,
    cv = "kfold"
  )

  expect_s4_class(fit, "GeneRankFit")
  expect_identical(fit@params$batch_col, "batch")
  expect_identical(fit@params$batch_covariates, "state")
})

test_that("rank_genes auto_confounds can switch cv to lobo and enable fold correction", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ranger")

  set.seed(6)

  expr <- matrix(rnorm(40 * 8), nrow = 40)
  rownames(expr) <- paste0("g", 1:40)
  colnames(expr) <- paste0("s", 1:8)

  # confounded on purpose: each batch mostly tied to one class
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = factor(c("A","A","A","B","B","B","B","A")),
      batch = factor(c("b1","b1","b1","b1","b2","b2","b2","b2"))
    )
  )

  fit <- suppressMessages(
    RFGeneRank::rank_genes(
      se,
      label_col = "state",
      batch_col = "batch",
      trees = 10,
      k = 2,
      n_top = 10,
      auto_confounds = TRUE,
      cv = "kfold"
    )
  )

  expect_s4_class(fit, "GeneRankFit")
  expect_true(fit@params$cv %in% c("kfold", "lobo"))
  expect_true(is.logical(fit@params$fold_batch_correction))
})

test_that("rank_genes returns zero PCA variance vector when fewer than two nonzero-variance genes remain", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ranger")

  expr <- matrix(1, nrow = 2, ncol = 8)
  rownames(expr) <- c("g1", "g2")
  colnames(expr) <- paste0("s", 1:8)
  expr[1, ] <- c(1,1,1,1,2,2,2,2)  # one variable gene only

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      state = factor(rep(c("A", "B"), each = 4))
    )
  )

  fit <- RFGeneRank::rank_genes(
    se,
    label_col = "state",
    trees = 10,
    k = 2,
    n_top = 0,
    auto_confounds = FALSE,
    cv = "kfold"
  )

  expect_s4_class(fit, "GeneRankFit")
  expect_equal(fit@var_pca, c(0, 0))
})