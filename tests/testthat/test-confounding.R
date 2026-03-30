test_that(".check_confounding returns severity and recommendations", {
  chk <- getFromNamespace(".check_confounding", "RFGeneRank")

  labels1 <- factor(rep(c("Control","Case"), each = 20))
  batch1  <- factor(rep(c("B1","B2"), times = 20))
  out1 <- chk(labels1, batch1, alpha = 0.01, v_thresh = 0.20)

  expect_true(is.list(out1))
  expect_true(all(c("p","cramerV","severity","recommendations") %in% names(out1)))
  expect_true(is.character(out1$severity))

  labels2 <- factor(c(rep("Control", 30), rep("Case", 30)))
  batch2  <- factor(c(rep("B1", 30), rep("B2", 30)))
  out2 <- chk(labels2, batch2, alpha = 0.01, v_thresh = 0.20)

  expect_true(is.list(out2))
  expect_true(out2$severity %in% c("moderate","high","severe","strong"))
  expect_true(is.list(out2$recommendations))
})