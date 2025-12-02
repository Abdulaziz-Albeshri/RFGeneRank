test_that("basic RF gene ranking workflow runs on toy data", {
  skip_on_cran()

  set.seed(1)
  expr <- matrix(
    rnorm(50 * 12),
    nrow = 50,
    dimnames = list(
      paste0("gene", seq_len(50)),
      paste0("sample", seq_len(12))
    )
  )

  state <- factor(rep(c("Control", "Case"), each = 6))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(expr = expr),
    colData = data.frame(
      state = state,
      row.names = colnames(expr)
    )
  )

  ranked <- rank_genes(
    se        = se,
    label_col = "state",
    trees     = 100,
    k         = 3,
    n_top     = 20,
    seed      = 1
  )

  # 1) Object exists
  expect_false(is.null(ranked))

  # 2) It has some content depending on its type
  if (is.data.frame(ranked)) {
    expect_gte(nrow(ranked), 1L)

  } else if (methods::is(ranked, "SummarizedExperiment")) {
    expect_gte(nrow(SummarizedExperiment::assay(ranked)), 1L)

  } else if (methods::is(ranked, "list")) {
    expect_gte(length(ranked), 1L)

  } else if (methods::is(ranked, "environment")) {
    expect_gte(length(ls(ranked)), 1L)

  } else {
    # Generic fallback for S4 or other objects:
    expect_gte(length(ranked), 1L)
  }
})
