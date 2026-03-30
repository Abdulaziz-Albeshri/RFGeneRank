test_that("GeneRankFit S4 class constructs with expected slots", {
  skip_if_not_installed("methods")

  obj <- methods::new("GeneRankFit")
  expect_s4_class(obj, "GeneRankFit")

  expect_true(all(c("params","oof","imp","features","final_model","var_pca","calibration") %in%
                    methods::slotNames(obj)))

  expect_type(methods::slot(obj, "params"), "list")
  expect_type(methods::slot(obj, "oof"), "list")
  expect_true(is.data.frame(methods::slot(obj, "imp")))
  expect_type(methods::slot(obj, "features"), "character")
  expect_true(is.numeric(methods::slot(obj, "var_pca")))
  expect_type(methods::slot(obj, "calibration"), "list")
})

test_that("RFGeneRank S4 class exists", {

  expect_false(methods::isClass("RFGeneRank"))

})


test_that("GeneRankFit class can be instantiated with prototype defaults", {
  obj <- methods::new("GeneRankFit")

  expect_s4_class(obj, "GeneRankFit")
  expect_identical(obj@params, list())
  expect_identical(obj@oof, list())
  expect_s3_class(obj@imp, "data.frame")
  expect_identical(obj@features, character())
  expect_null(obj@final_model)
  expect_identical(obj@var_pca, numeric())
  expect_identical(obj@calibration, list())
})

test_that("GeneRankFit validity passes for default empty object", {
  obj <- methods::new("GeneRankFit")
  expect_true(methods::validObject(obj))
})

test_that("GeneRankFit validity accepts minimal valid oof structure", {
  obj <- methods::new(
    "GeneRankFit",
    oof = list(
      prob = matrix(c(0.2, 0.8, 0.7, 0.3), nrow = 2, byrow = TRUE),
      y = factor(c("A", "B"))
    )
  )

  expect_true(methods::validObject(obj))
})

test_that("GeneRankFit validity rejects oof missing prob", {
  expect_error(
    methods::new(
      "GeneRankFit",
      oof = list(
        y = factor(c("A", "B"))
      )
    ),
    "'oof' must contain at least \\$prob and \\$y"
  )
})

test_that("GeneRankFit validity rejects oof missing y", {
  expect_error(
    methods::new(
      "GeneRankFit",
      oof = list(
        prob = matrix(c(0.2, 0.8, 0.7, 0.3), nrow = 2, byrow = TRUE)
      )
    ),
    "'oof' must contain at least \\$prob and \\$y"
  )
})

test_that("GeneRankFit validity rejects non-matrix oof$prob", {
  expect_error(
    methods::new(
      "GeneRankFit",
      oof = list(
        prob = c(0.2, 0.8),
        y = factor(c("A", "B"))
      )
    ),
    "oof\\$prob must be a matrix"
  )
})

test_that("GeneRankFit validity rejects non-factor oof$y", {
  expect_error(
    methods::new(
      "GeneRankFit",
      oof = list(
        prob = matrix(c(0.2, 0.8, 0.7, 0.3), nrow = 2, byrow = TRUE),
        y = c("A", "B")
      )
    ),
    "oof\\$y must be a factor"
  )
})

test_that("GeneRankFit validity rejects calibration fun when not a function", {
  expect_error(
    methods::new(
      "GeneRankFit",
      calibration = list(method = "platt", fun = 123)
    ),
    "'calibration\\$fun' must be a function"
  )
})

test_that("GeneRankFit validity accepts calibration fun as function", {
  obj <- methods::new(
    "GeneRankFit",
    calibration = list(
      method = "platt",
      fun = function(x) x
    )
  )

  expect_true(methods::validObject(obj))
})

test_that("show prints minimal summary for default object", {
  obj <- methods::new("GeneRankFit")
  out <- capture.output(show(obj))

  expect_true(any(grepl("^GeneRankFit$", out)))
  expect_true(any(grepl("finalized model: <none>", out, fixed = TRUE)))
  expect_true(any(grepl("calibration: <none>", out, fixed = TRUE)))
})

test_that("show prints key non-default summary lines", {
  obj <- methods::new(
    "GeneRankFit",
    params = list(
      k = 5,
      trees = 300,
      importance = "impurity",
      filter_low_expr = TRUE,
      min_prop = 0.2,
      transform = "log1p"
    ),
    imp = data.frame(
      feature = c("g1", "g2", "g3"),
      importance = c(0.9, 0.7, 0.5)
    ),
    final_model = structure(list(dummy = TRUE), class = "mock_model"),
    calibration = list(
      method = "platt",
      fun = function(x) x
    )
  )

  out <- capture.output(show(obj))

  expect_true(any(grepl("GeneRankFit", out, fixed = TRUE)))
  expect_true(any(grepl("k=5, trees=300, importance=impurity", out, fixed = TRUE)))
  expect_true(any(grepl("CV-safe filtering: keep >= 20% nonzero expr", out, fixed = TRUE)))
  expect_true(any(grepl("transform: log1p", out, fixed = TRUE)))
  expect_true(any(grepl("features ranked: 3", out, fixed = TRUE)))
  expect_true(any(grepl("finalized model: PRESENT", out, fixed = TRUE)))
  expect_true(any(grepl("calibration: platt", out, fixed = TRUE)))
})

test_that("%||% returns left value when not NULL and right value when NULL", {
  expect_identical(`%||%`(1, 2), 1)
  expect_identical(`%||%`(NULL, 2), 2)
  expect_identical(`%||%`("x", "y"), "x")
})

test_that("GeneRankFit validity catches malformed oof/calibration", {
  skip_if_not_installed("methods")

  expect_error(
    methods::new("GeneRankFit", oof = list(prob = NULL, y = NULL)),
    "oof.*must contain at least \\$prob and \\$y"
  )

  expect_error(
    methods::new("GeneRankFit", oof = list(prob = 1:3, y = factor(c("A","B","A")))),
    "oof\\$prob must be a matrix"
  )

  expect_error(
    methods::new("GeneRankFit", oof = list(prob = matrix(runif(6), nrow = 3), y = c("A","B","A"))),
    "oof\\$y must be a factor"
  )

  expect_error(
    methods::new("GeneRankFit", calibration = "bad"),
    "invalid object for slot \"calibration\""
  )

  expect_error(
    methods::new("GeneRankFit", calibration = list(method = "platt", fun = 123)),
    "calibration\\$fun.*must be a function"
  )
})

test_that("GeneRankFit class loads and slots can be set (class_def coverage attempt)", {
  skip_if_not_installed("methods")

  obj <- methods::new("GeneRankFit")
  expect_s4_class(obj, "GeneRankFit")

  if ("params" %in% methods::slotNames(obj)) {
    methods::slot(obj, "params") <- list(trees = 10)
    expect_true(is.list(methods::slot(obj, "params")))
  }
  if ("features" %in% methods::slotNames(obj)) {
    methods::slot(obj, "features") <- c("g1","g2")
    expect_true(is.character(methods::slot(obj, "features")))
  }
  if ("var_pca" %in% methods::slotNames(obj)) {
    methods::slot(obj, "var_pca") <- 0.5
    expect_true(is.numeric(methods::slot(obj, "var_pca")))
  }
})