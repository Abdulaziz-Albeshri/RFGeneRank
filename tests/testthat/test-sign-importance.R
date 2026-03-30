test_that("sign_importance(method='mean') works on small matrices", {
  set.seed(1)
  X <- matrix(abs(rnorm(6 * 10)) + 0.01, nrow = 6)
  rownames(X) <- paste0("s", 1:6)
  colnames(X) <- paste0("g", 1:10)
  y <- factor(rep(c("Control","Case"), each = 3))

  fit <- list(importance = setNames(runif(10), colnames(X)))

  tab <- sign_importance(fit = fit, X = X, y = y, method = "mean")
  expect_true(is.data.frame(tab))
  expect_true(all(c("gene","importance") %in% colnames(tab)))
})

test_that("sign_importance errors when y is unavailable", {
  X <- matrix(rnorm(4 * 5), nrow = 4)
  colnames(X) <- paste0("g", 1:5)

  fit <- list(
    importance = setNames(c(5, 4, 3, 2, 1), colnames(X))
  )

  expect_error(
    sign_importance(fit = fit, X = X, y = NULL, method = "mean"),
    regexp = "required if it cannot be extracted|Dimensions of X and y do not align"
  )
})

test_that("sign_importance errors when X and y dimensions do not align explicitly", {
  X <- matrix(rnorm(4 * 5), nrow = 4)
  colnames(X) <- paste0("g", 1:5)
  y <- factor(c("A", "A", "B"))

  fit <- list(
    importance = setNames(c(5, 4, 3, 2, 1), colnames(X))
  )

  expect_error(
    sign_importance(fit = fit, X = X, y = y, method = "mean"),
    "Dimensions of X and y do not align"
  )
})

test_that("sign_importance returns structured output with method mean", {
  X <- matrix(rnorm(4 * 5), nrow = 4)
  colnames(X) <- paste0("g", 1:5)
  y <- factor(c("Control", "Control", "Case", "Case"))

  fit <- list(
    importance = setNames(c(0.8, 0.6, 0.4, 0.2, 0.1), colnames(X))
  )

  res <- suppressWarnings(
    sign_importance(
      fit = fit,
      X = X,
      y = y,
      method = "mean"
    )
  )

  expect_true(is.data.frame(res))
  expect_true(all(c(
    "gene", "importance", "direction", "signed_importance",
    "mean_case", "mean_ctrl", "mean_diff", "log2FC", "shap_dir"
  ) %in% names(res)))
  expect_equal(nrow(res), 5)
})

test_that("sign_importance works with method de", {
  X <- matrix(abs(rnorm(4 * 5)), nrow = 4)
  colnames(X) <- paste0("g", 1:5)
  y <- factor(c("Control", "Control", "Case", "Case"))

  fit <- list(
    importance = setNames(c(0.8, 0.6, 0.4, 0.2, 0.1), colnames(X))
  )

  de_table <- data.frame(
    gene = colnames(X),
    log2FC = c(1, -1, 0.5, -0.5, 0)
  )

  res <- sign_importance(
    fit = fit,
    X = X,
    y = y,
    method = "de",
    de_table = de_table
  )

  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 5)
  expect_true(all(res$gene %in% colnames(X)))
})

test_that("sign_importance errors for de method without proper de_table", {
  X <- matrix(abs(rnorm(4 * 5)), nrow = 4)
  colnames(X) <- paste0("g", 1:5)
  y <- factor(c("Control", "Control", "Case", "Case"))

  fit <- list(
    importance = setNames(c(0.8, 0.6, 0.4, 0.2, 0.1), colnames(X))
  )

  expect_error(
    sign_importance(
      fit = fit,
      X = X,
      y = y,
      method = "de",
      de_table = data.frame(gene = colnames(X))
    ),
    "For method='de', provide de_table with columns: gene, log2FC"
  )
})

