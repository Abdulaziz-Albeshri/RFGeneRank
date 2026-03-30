test_that("calibration getter/setter and apply_calibration are stable", {
  fit <- methods::new("GeneRankFit")

  p <- c(0.1, 0.4, 0.9)
  expect_equal(apply_calibration(fit, p), p)

  calibration(fit) <- list(method = "custom", fun = function(x) x^0.8)
  out <- apply_calibration(fit, p)

  expect_true(is.numeric(out))
  expect_length(out, length(p))
  expect_true(all(out >= 0 & out <= 1))
})

test_that("apply_calibration returns original probabilities when unset", {
  fit <- methods::new("GeneRankFit")
  p <- c(0.1, 0.5, 0.9)
  expect_equal(apply_calibration(fit, p), p)
})

test_that("calibrate_oof(platt) stores fun and changes output", {
  set.seed(1)
  n <- 30
  y <- factor(rep(c("Control","Case"), each = n/2))

  p_case <- c(runif(n/2, 0.05, 0.35), runif(n/2, 0.65, 0.95))
  prob <- cbind(Control = 1 - p_case, Case = p_case)
  rownames(prob) <- paste0("s", 1:n)
  names(y) <- rownames(prob)

  fit <- methods::new("GeneRankFit")
  slot(fit, "oof") <- list(prob = prob, y = y)

  fit2 <- calibrate_oof(fit, method = "platt")
  cal <- calibration(fit2)

  expect_true(is.list(cal))
  expect_true(is.function(cal$fun))
  expect_identical(cal$method, "platt")

  x <- c(0.1, 0.3, 0.7, 0.9)
  expect_false(isTRUE(all.equal(apply_calibration(fit2, x), x)))
})

test_that("calibrate_oof(isotonic) returns bounded outputs", {
  set.seed(2)
  n <- 40
  y <- factor(rep(c("Control","Case"), length.out = n))

  p_case <- pmin(pmax(rnorm(n, mean = ifelse(y=="Case", 0.75, 0.25), sd = 0.12), 0), 1)
  prob <- cbind(Control = 1 - p_case, Case = p_case)
  rownames(prob) <- paste0("s", 1:n)
  names(y) <- rownames(prob)

  fit <- methods::new("GeneRankFit")
  slot(fit, "oof") <- list(prob = prob, y = y)

  fit2 <- calibrate_oof(fit, method = "isotonic")

  z <- apply_calibration(fit2, seq(0, 1, length.out = 9))
  expect_true(all(z >= 0 & z <= 1))
})

test_that("calibrate_oof errors on missing rownames and handles single-class input", {
  # missing rownames -> should error
  prob <- matrix(runif(20), nrow = 10)
  colnames(prob) <- c("Control", "Case")
  y <- factor(rep(c("Control", "Case"), length.out = 10))

  fit <- methods::new("GeneRankFit")
  slot(fit, "oof") <- list(prob = prob, y = y)

  expect_error(
    calibrate_oof(fit, "platt"),
    "rownames"
  )

  # single class
  prob2 <- cbind(Control = rep(0.9, 10), Case = rep(0.1, 10))
  rownames(prob2) <- paste0("s", 1:10)
  y2 <- factor(rep("Control", 10), levels = c("Control", "Case"))
  names(y2) <- rownames(prob2)

  fit2 <- methods::new("GeneRankFit")
  slot(fit2, "oof") <- list(prob = prob2, y = y2)

  # check the real current behaviour
  expect_error(
    calibrate_oof(fit2, "isotonic"),
    "both classes|two classes|binary",
    ignore.case = TRUE
  )
})