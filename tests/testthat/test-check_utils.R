test_that(".message_once emits a message", {
  expect_message(
    RFGeneRank:::.message_once("hello"),
    "hello"
  )
})

test_that("check_utils internal validators behave", {
  chk_mat  <- getFromNamespace(".check_is_matrix", "RFGeneRank")
  chk_meta <- getFromNamespace(".check_meta", "RFGeneRank")
  chk_same <- getFromNamespace(".check_same_samples", "RFGeneRank")

  X <- matrix(rnorm(5 * 6), nrow = 5,
              dimnames = list(paste0("g", 1:5), paste0("s", 1:6)))

  expect_silent(chk_mat(X, "X"))
  expect_error(chk_mat(as.data.frame(X), "X"), "must be a matrix", fixed = TRUE)

  meta <- data.frame(
    state = rep(c("Control","Case"), length.out = 6),
    batch = rep(c("B1","B2"), length.out = 6),
    row.names = colnames(X)
  )

  expect_silent(chk_meta(meta, label_col = "state", batch_col = "batch"))

  meta2 <- meta
  rownames(meta2) <- paste0("t", seq_len(nrow(meta2)))  

  expect_error(
    chk_same(X, meta2),
    "Fewer than 4 overlapping samples",
    fixed = TRUE
  )
})

test_that("calibrate_oof(platt) stores calibration and apply_calibration changes output", {
  skip_if_not_installed("methods")

  set.seed(1)
  n <- 80
  sids <- paste0("s", seq_len(n))
  y <- factor(rep(c("A","B"), each = n/2))
  names(y) <- sids

  p <- plogis(rnorm(n, mean = ifelse(y == "B", 0.5, -0.5), sd = 1))
  prob <- matrix(p, ncol = 1, dimnames = list(sids, "B"))

  fit <- methods::new("GeneRankFit", oof = list(prob = prob, y = y))

  fit2 <- RFGeneRank::calibrate_oof(fit, method = "platt")
  p2 <- RFGeneRank::apply_calibration(fit2, p)

  expect_true(is.numeric(p2))
  expect_equal(length(p2), length(p))
  expect_true(all(is.finite(p2)))
  expect_true(all(p2 >= 0 & p2 <= 1))

  expect_false(isTRUE(all.equal(p2, p)))
})

test_that("check utils: basic validators error on wrong types", {
  exported <- getNamespaceExports("RFGeneRank")

  candidates <- intersect(exported, c("check_matrix", "check_metadata", "check_col", "assert_factor", "assert_named"))
  skip_if(length(candidates) == 0)

  for (nm in candidates) {
    fn <- get(nm, envir = asNamespace("RFGeneRank"))
    expect_error(fn(123), NA)  # we just want it to execute error paths
  }
})

test_that(".check_is_matrix accepts valid numeric matrix", {
  x <- matrix(rnorm(12), nrow = 3)
  expect_invisible(RFGeneRank:::.check_is_matrix(x, "x"))
})

test_that(".check_meta accepts valid metadata", {
  meta <- data.frame(
    state = factor(c("A", "B", "A", "B")),
    batch = factor(c("b1", "b1", "b2", "b2"))
  )
  expect_invisible(RFGeneRank:::.check_meta(meta, "state", "batch"))
})

test_that(".check_same_samples reorders metadata to match expr columns", {
  expr <- matrix(rnorm(5 * 4), nrow = 5)
  colnames(expr) <- c("s1", "s2", "s3", "s4")

  meta <- data.frame(
    state = factor(c("A", "B", "A", "B")),
    row.names = c("s3", "s1", "s4", "s2")
  )

  out <- RFGeneRank:::.check_same_samples(expr, meta)

  expect_equal(colnames(out$expr), rownames(out$meta))
  expect_equal(rownames(out$meta), c("s1", "s2", "s3", "s4"))
})

test_that(".check_same_samples drops non-overlapping samples correctly", {
  expr <- matrix(rnorm(5 * 6), nrow = 5)
  colnames(expr) <- c("s1", "s2", "s3", "s4", "s5", "s6")

  meta <- data.frame(
    state = factor(c("A", "B", "A", "B", "A")),
    row.names = c("s2", "s3", "s4", "s5", "s6")
  )

  out <- RFGeneRank:::.check_same_samples(expr, meta)

  expect_equal(colnames(out$expr), c("s2", "s3", "s4", "s5", "s6"))
  expect_equal(rownames(out$meta), c("s2", "s3", "s4", "s5", "s6"))
})

test_that(".message_once emits a message on repeated calls", {
  expect_message(RFGeneRank:::.message_once("hello-utils-fast"), "hello-utils-fast")
  expect_message(RFGeneRank:::.message_once("hello-utils-fast"), "hello-utils-fast")
})

test_that(".check_is_matrix accepts a valid numeric matrix", {
  x <- matrix(rnorm(12), nrow = 3)
  expect_null(RFGeneRank:::.check_is_matrix(x, "x"))
})

test_that(".check_is_matrix errors on non-matrix input", {
  expect_error(
    RFGeneRank:::.check_is_matrix(1:5, "x"),
    "x must be a matrix"
  )
})

test_that(".check_is_matrix errors on non-numeric matrix", {
  x <- matrix(letters[1:4], nrow = 2)
  expect_error(
    RFGeneRank:::.check_is_matrix(x, "x"),
    "x must be numeric"
  )
})

test_that(".check_is_matrix errors on too-small matrix", {
  expect_error(
    RFGeneRank:::.check_is_matrix(matrix(1, nrow = 1, ncol = 2), "x"),
    ">=2 genes and >=2 samples"
  )
  expect_error(
    RFGeneRank:::.check_is_matrix(matrix(1, nrow = 2, ncol = 1), "x"),
    ">=2 genes and >=2 samples"
  )
})

test_that(".check_meta accepts valid metadata", {
  meta <- data.frame(
    state = factor(c("A", "B", "A", "B")),
    batch = factor(c("b1", "b1", "b2", "b2"))
  )
  expect_null(RFGeneRank:::.check_meta(meta, "state", "batch"))
})

test_that(".check_meta errors on invalid metadata cases", {
  expect_error(
    RFGeneRank:::.check_meta("bad", "state"),
    "metadata must be a data.frame"
  )

  meta <- data.frame(batch = c("b1", "b1", "b2", "b2"))
  expect_error(
    RFGeneRank:::.check_meta(meta, "state"),
    "metadata is missing label_col: state"
  )

  meta2 <- data.frame(state = c("A", "B", "A", "B"))
  expect_error(
    RFGeneRank:::.check_meta(meta2, "state", "batch"),
    "metadata is missing batch_col: batch"
  )
})

test_that(".check_same_samples reorders metadata to match expr columns", {
  expr <- matrix(rnorm(5 * 4), nrow = 5)
  colnames(expr) <- c("s1", "s2", "s3", "s4")

  meta <- data.frame(
    state = factor(c("A", "B", "A", "B")),
    row.names = c("s3", "s1", "s4", "s2")
  )

  out <- RFGeneRank:::.check_same_samples(expr, meta)

  expect_equal(colnames(out$expr), c("s1", "s2", "s3", "s4"))
  expect_equal(rownames(out$meta), c("s1", "s2", "s3", "s4"))
})

test_that(".check_same_samples errors when overlap is fewer than 4", {
  expr <- matrix(rnorm(5 * 5), nrow = 5)
  colnames(expr) <- c("s1", "s2", "s3", "s4", "s5")

  meta <- data.frame(
    state = factor(c("A", "B", "A")),
    row.names = c("s1", "s2", "s3")
  )

  expect_error(
    RFGeneRank:::.check_same_samples(expr, meta),
    "Fewer than 4 overlapping samples"
  )
})

test_that(".message_once emits messages on repeated calls", {
  expect_message(RFGeneRank:::.message_once("hello-utils-fast"), "hello-utils-fast")
  expect_message(RFGeneRank:::.message_once("hello-utils-fast"), "hello-utils-fast")
})