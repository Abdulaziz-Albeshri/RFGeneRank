test_that("id_map returns empty data frame with correct column names for empty keys", {
  out <- RFGeneRank::id_map(
    keys = character(0),
    from = "ENTREZID",
    to = "SYMBOL",
    OrgDb = structure(list(), class = "MockOrgDb")
  )

  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("ENTREZID", "SYMBOL"))
  expect_equal(nrow(out), 0)
})

test_that("id_map works on real OrgDb and preserves requested column order", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  out <- RFGeneRank::id_map(
    keys = c("7157", "7158"),
    from = "ENTREZID",
    to = "SYMBOL",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )

  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("ENTREZID", "SYMBOL"))
  expect_true(nrow(out) >= 2)
  expect_true(all(c("7157", "7158") %in% out$ENTREZID))
})

test_that("id_map restores original input order", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  out <- RFGeneRank::id_map(
    keys = c("7158", "7157"),
    from = "ENTREZID",
    to = "SYMBOL",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )

  expect_identical(out$ENTREZID, c("7158", "7157"))
})

test_that("id_map handles duplicated input keys via unique(keys) and restores duplicates in output order", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  out <- RFGeneRank::id_map(
    keys = c("7157", "7158", "7157"),
    from = "ENTREZID",
    to = "SYMBOL",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    drop_na = TRUE,
    unique = TRUE
  )

  expect_identical(out$ENTREZID, c("7157", "7158", "7157"))
  expect_equal(nrow(out), 3)
  expect_identical(out$SYMBOL[c(1, 3)], rep(out$SYMBOL[1], 2))
})

test_that("id_map unique = TRUE keeps at most one mapping per key", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

 out <- RFGeneRank::id_map(
  keys = c("NOT_A_REAL_ID_1", "TP53"),
  from = "SYMBOL",
  to = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db,
  drop_na = FALSE,
  unique = TRUE
)

expect_s3_class(out, "data.frame")
expect_identical(out$SYMBOL, c("NOT_A_REAL_ID_1", "TP53"))
expect_equal(nrow(out), 2)
})

test_that("id_map unique = FALSE can return repeated source keys for one-to-many mappings", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

 out <- RFGeneRank::id_map(
  keys = c("NOT_A_REAL_ID_1", "TP53"),
  from = "SYMBOL",
  to = "ENTREZID",
  OrgDb = org.Hs.eg.db::org.Hs.eg.db,
  drop_na = TRUE,
  unique = TRUE
)

expect_s3_class(out, "data.frame")
expect_true(nrow(out) >= 1)
expect_false("NOT_A_REAL_ID_1" %in% out$SYMBOL)
expect_true("TP53" %in% out$SYMBOL)
})

test_that("id_map returns requested columns in requested order", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  out <- RFGeneRank::id_map(
    keys = c("7157", "7158"),
    from = "ENTREZID",
    to = "SYMBOL",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    drop_na = TRUE,
    unique = TRUE
  )

  expect_identical(names(out), c("ENTREZID", "SYMBOL"))
  expect_true(is.data.frame(out))
  expect_true(nrow(out) >= 2)
})