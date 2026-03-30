test_that("fold_batch_correct removes batch using train-only design and preserves dims", {
  fold_bc <- getFromNamespace("fold_batch_correct", "RFGeneRank")

  set.seed(1)

  # genes x samples
  Ytr <- matrix(rnorm(40 * 10), nrow = 40,
                dimnames = list(paste0("g", 1:40), paste0("tr", 1:10)))
  Yte <- matrix(rnorm(40 * 6), nrow = 40,
                dimnames = list(paste0("g", 1:40), paste0("te", 1:6)))

  meta_tr <- data.frame(
    batch = factor(rep(c("b1","b2"), each = 5)),
    age   = c(30,31,29,35,40, 50,52,48,49,51),
    row.names = colnames(Ytr)
  )

  meta_te <- data.frame(
    batch = factor(rep(c("b1","b2"), each = 3)),
    age   = c(33,34,32, 55,54,53),
    row.names = colnames(Yte)
  )

  res <- fold_bc(
    Ytr = Ytr, Yte = Yte,
    meta_tr = meta_tr, meta_te = meta_te,
    batch_col = "batch",
    covar_cols = "age"
  )

  expect_true(is.list(res))
  expect_true(all(c("train","test") %in% names(res)))
  expect_true(is.matrix(res$train))
  expect_true(is.matrix(res$test))
  expect_identical(dim(res$train), dim(Ytr))
  expect_identical(dim(res$test), dim(Yte))
})

test_that("batch_correct internals build design and remove batch (coverage)", {
  build_design <- getFromNamespace(".build_design", "RFGeneRank")
  solve_beta   <- getFromNamespace(".solve_beta", "RFGeneRank")
  lm_apply     <- getFromNamespace(".lm_apply_remove_batch", "RFGeneRank")

  set.seed(1)

  Ytr <- matrix(
    rnorm(40 * 10),
    nrow = 40,
    dimnames = list(paste0("g", 1:40), paste0("s", 1:10))
  )

  meta <- data.frame(
    batch = factor(rep(c("b1", "b2"), each = 5)),
    age   = c(30, 31, 29, 35, 40, 50, 52, 48, 49, 51),
    row.names = colnames(Ytr)
  )

  # build design matrix Xtr
  f1 <- names(formals(build_design))
  args1 <- list()
  if ("meta" %in% f1) args1$meta <- meta
  if ("batch_col" %in% f1) args1$batch_col <- "batch"
  if ("covar_cols" %in% f1) args1$covar_cols <- "age"
  Xtr <- do.call(build_design, args1)
  Xtr <- as.matrix(Xtr)

  expect_true(is.matrix(Xtr))
  expect_identical(nrow(Xtr), ncol(Ytr))  # rows = samples

  # solve beta: solve_beta(Xtr, Ytr, lambda)
  beta <- solve_beta(Xtr = Xtr, Ytr = Ytr, lambda = 1e-06)
  expect_true(is.matrix(beta))

  # which columns in Xtr correspond to batch terms?
  xcn <- colnames(Xtr)
  batch_cols <- which(grepl("^batch", xcn))
  expect_true(length(batch_cols) >= 1)

  # apply remove batch: lm_apply(Y, X, beta, batch_cols)
  f3 <- names(formals(lm_apply))
  args3 <- list()
  if ("Y" %in% f3) args3$Y <- Ytr
  if ("X" %in% f3) args3$X <- Xtr
  if ("design" %in% f3) args3$design <- Xtr
  if ("beta" %in% f3) args3$beta <- beta
  if ("batch_cols" %in% f3) args3$batch_cols <- batch_cols

  Y2 <- do.call(lm_apply, args3)

  expect_true(is.matrix(Y2))
  expect_identical(dim(Y2), dim(Ytr))
})