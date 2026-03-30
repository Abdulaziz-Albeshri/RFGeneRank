test_that(".compute_class_weights errors when y is NULL", {
  expect_error(
    RFGeneRank:::.compute_class_weights(NULL),
    "`y` is NULL"
  )
})

test_that(".compute_class_weights errors when fewer than two classes remain", {
  expect_error(
    RFGeneRank:::.compute_class_weights(factor(rep("A", 5))),
    "needs at least two classes"
  )

  expect_error(
    RFGeneRank:::.compute_class_weights(c("A", "A", "A")),
    "needs at least two classes"
  )
})

test_that(".compute_class_weights removes NAs when na.rm = TRUE", {
  y <- c("A", "A", NA, "B", "B", NA)

  w <- RFGeneRank:::.compute_class_weights(y, na.rm = TRUE)

  expect_true(is.numeric(w))
  expect_true(all(c("A", "B") %in% names(w)))
  expect_equal(length(w), 2)
  expect_equal(unname(mean(w)), 1, tolerance = 1e-8)
})

test_that(".compute_class_weights keeps balanced classes at weight 1", {
  y <- factor(c("A", "A", "B", "B"))

  w <- RFGeneRank:::.compute_class_weights(y)

  expect_true(is.numeric(w))
  expect_equal(length(w), 2)
  expect_equal(unname(w["A"]), 1, tolerance = 1e-8)
  expect_equal(unname(w["B"]), 1, tolerance = 1e-8)
  expect_equal(unname(mean(w)), 1, tolerance = 1e-8)
})

test_that(".compute_class_weights gives larger weight to rarer class", {
  y <- factor(c(rep("A", 6), rep("B", 2)))

  w <- RFGeneRank:::.compute_class_weights(y)

  expect_true(is.numeric(w))
  expect_true(w["B"] > w["A"])
  expect_equal(unname(mean(w)), 1, tolerance = 1e-8)
})

test_that(".compute_class_weights coerces non-factor input and returns named vector", {
  y <- c("case", "case", "case", "ctrl", "ctrl")

  w <- RFGeneRank:::.compute_class_weights(y)

  expect_true(is.numeric(w))
  expect_true(!is.null(names(w)))
  expect_setequal(names(w), c("case", "ctrl"))
  expect_equal(unname(mean(w)), 1, tolerance = 1e-8)
})

test_that(".compute_class_weights respects na.rm = FALSE when NAs are present", {
  y <- factor(c("A", "A", NA, "B", "B"), exclude = NULL)

  w <- RFGeneRank:::.compute_class_weights(y, na.rm = FALSE)

  expect_true(is.numeric(w))
  expect_true(length(w) >= 2)
  expect_equal(unname(mean(w)), 1, tolerance = 1e-8)
})

test_that("class weights returns named numeric", {
  y <- factor(c(rep("A", 8), rep("B", 2)))

  fn <- getFromNamespace(".compute_class_weights","RFGeneRank")
  w <- fn(y)

  expect_true(is.numeric(w))
  expect_true(!is.null(names(w)))
  expect_true(all(levels(y) %in% names(w)))
})

test_that("ComBat train/apply helpers preserve dims", {
  skip_if_not_installed("sva")
  # genes x samples
  Xtr <- matrix(rnorm(50), nrow = 10,
                dimnames = list(paste0("g",1:10), paste0("tr",1:5)))
  Xte <- matrix(rnorm(30), nrow = 10,
                dimnames = list(paste0("g",1:10), paste0("te",1:3)))

  batch_tr <- factor(c("b1","b1","b2","b2","b2"))
  mod_tr <- model.matrix(~1, data = data.frame(batch_tr))

  combat_fit <- getFromNamespace("rfgr_combat_fit_train","RFGeneRank")
  combat_apply <- getFromNamespace("rfgr_combat_apply_test","RFGeneRank")

  fit <- combat_fit(X_tr = Xtr, batch_tr = batch_tr, mod_tr = mod_tr)

  Xtr2 <- combat_apply(Xtr, batch_te = batch_tr, mod_te = mod_tr, estimates = fit$estimates)
  batch_te <- factor(c("b1","b2","b2"))
  mod_te <- model.matrix(~1, data = data.frame(batch_te))
  Xte2 <- combat_apply(Xte, batch_te = batch_te, mod_te = mod_te, estimates = fit$estimates)

  expect_identical(dim(Xtr2), dim(Xtr))
  expect_identical(dim(Xte2), dim(Xte))
})
