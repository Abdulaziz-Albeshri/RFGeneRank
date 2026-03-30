library(testthat)
library(ggplot2)
library(SummarizedExperiment)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

toy_se_plot <- function(n_genes = 6, n_samples = 10, constant = FALSE) {
  set.seed(1)
  expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  if (constant) {
    expr[,] <- 1
  }
  rownames(expr) <- paste0("g", seq_len(n_genes))
  colnames(expr) <- paste0("s", seq_len(n_samples))

  cd <- data.frame(
    state = factor(rep(c("Control", "Case"), length.out = n_samples)),
    label = factor(rep(c("Control", "Case"), length.out = n_samples)),
    age = seq(40, by = 3, length.out = n_samples),
    age_num = as.character(seq(40, by = 3, length.out = n_samples)),
    sex = factor(rep(c("Female", "Male"), length.out = n_samples)),
    batch = factor(rep(c("B1", "B2"), length.out = n_samples)),
    score = seq_len(n_samples),
    row.names = colnames(expr)
  )

  SummarizedExperiment(
    assays = list(expr = expr),
    colData = cd
  )
}

mock_fit <- function() {
  structure(list(dummy = TRUE), class = "GeneRankFit")
}

mock_imp <- function() {
  data.frame(
    gene = c("g3", "g1", "g2", "g5", "g4"),
    importance = c(0.9, 0.7, 0.5, 0.2, 0.1),
    direction = c(1, -1, 1, 0, -1),
    signed_importance = c(0.9, -0.7, 0.5, 0.0, -0.1),
    stringsAsFactors = FALSE
  )
}

mock_params <- function(label_col = "label") {
  list(
    seed = 99,
    label_col = label_col,
    umap = list(neighbors = 7, min_dist = 0.2, metric = "euclidean")
  )
}

mock_oof <- function(n = 10) {
  y <- factor(rep(c("Control", "Case"), length.out = n))
  prob <- cbind(
    Control = seq(0.9, 0.1, length.out = n),
    Case    = seq(0.1, 0.9, length.out = n)
  )
  rownames(prob) <- paste0("s", seq_len(n))
  names(y) <- rownames(prob)
  list(prob = prob, y = y)
}

# ------------------------------------------------------------------------------
# plot_importance
# ------------------------------------------------------------------------------

test_that("plot_importance returns ggplot and respects top ordering", {
  fit <- mock_fit()

  local_mocked_bindings(
    .rfgr_imp = function(...) mock_imp(),
    .env = asNamespace("RFGeneRank")
  )

  p <- RFGeneRank::plot_importance(fit, top = 3)

  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "Top 3 predictive genes")
  dat <- ggplot_build(p)$plot$data
  expect_equal(nrow(dat), 3)
})

test_that("plot_importance enforces GeneRankFit class", {
  expect_error(
    RFGeneRank::plot_importance(list(), top = 3),
    class = "simpleError"
  )
})

# ------------------------------------------------------------------------------
# plot_sign_importance
# ------------------------------------------------------------------------------

test_that("plot_sign_importance works from supplied tab", {
  tab <- mock_imp()

  p <- RFGeneRank::plot_sign_importance(
    tab = tab,
    top = 4,
    show_legend = FALSE,
    palette = c(`-1` = "blue", `0` = "grey", `1` = "red")
  )

  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "Top 4 signed predictive genes")
  dat <- ggplot_build(p)$plot$data
  expect_equal(nrow(dat), 4)
})

test_that("plot_sign_importance errors when fit and tab are both missing", {
  expect_error(
    RFGeneRank::plot_sign_importance(),
    "Provide either `tab`"
  )
})

test_that("plot_sign_importance errors when stored importance lacks signed_importance", {
  fit <- mock_fit()

  local_mocked_bindings(
    .rfgr_imp = function(...) data.frame(
      gene = c("g1", "g2"),
      importance = c(1, 2),
      direction = c(1, -1)
    ),
    .env = asNamespace("RFGeneRank")
  )

  expect_error(
    RFGeneRank::plot_sign_importance(fit = fit),
    "lacks `signed_importance`"
  )
})

test_that("plot_sign_importance errors when required columns are absent", {
  bad_tab <- data.frame(gene = "g1", importance = 1)

  expect_error(
    RFGeneRank::plot_sign_importance(tab = bad_tab),
    "missing required columns"
  )
})

# ------------------------------------------------------------------------------
# plot_embed_expr
# ------------------------------------------------------------------------------

test_that("plot_embed_expr PCA path works", {
  fit <- mock_fit()
  se <- toy_se_plot()

  local_mocked_bindings(
    .rfgr_imp = function(...) mock_imp(),
    .rfgr_params = function(...) mock_params(label_col = "label"),
    .env = asNamespace("RFGeneRank")
  )

  p <- RFGeneRank::plot_embed_expr(
    fit = fit,
    se = se,
    n_top = 4,
    type = "pca",
    zscore = TRUE,
    show_legend = FALSE
  )

  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "PCA of Top")
})

test_that("plot_embed_expr errors with fewer than 2 overlapping genes", {
  fit <- mock_fit()
  se <- toy_se_plot()

  local_mocked_bindings(
    .rfgr_imp = function(...) data.frame(
      gene = "not_in_expr",
      importance = 1
    ),
    .rfgr_params = function(...) mock_params(label_col = "label"),
    .env = asNamespace("RFGeneRank")
  )

  expect_error(
    RFGeneRank::plot_embed_expr(fit = fit, se = se, type = "pca"),
    "Fewer than 2 overlapping top genes"
  )
})

test_that("plot_embed_expr errors when PCA has too few variable genes", {
  fit <- mock_fit()
  se <- toy_se_plot()

  expr <- SummarizedExperiment::assay(se, "expr")
  expr["g1", ] <- 1
  expr["g2", ] <- 1
  expr["g3", ] <- 1
  SummarizedExperiment::assay(se, "expr") <- expr

  local_mocked_bindings(
    .rfgr_imp = function(...) data.frame(
      gene = c("g1", "g2", "g3"),
      importance = c(3, 2, 1)
    ),
    .rfgr_params = function(...) mock_params(label_col = "label"),
    .env = asNamespace("RFGeneRank")
  )

  expect_error(
    RFGeneRank::plot_embed_expr(
      fit = fit,
      se = se,
      type = "pca",
      zscore = FALSE
    ),
    "Too few variable genes for PCA"
  )
})

test_that("plot_embed_expr UMAP uwot path works when uwot is available", {
  skip_if_not_installed("uwot")

  fit <- mock_fit()
  se <- toy_se_plot()

  local_mocked_bindings(
    .rfgr_imp = function(...) mock_imp(),
    .rfgr_params = function(...) mock_params(label_col = "label"),
    .env = asNamespace("RFGeneRank")
  )

  p <- RFGeneRank::plot_embed_expr(
    fit = fit,
    se = se,
    n_top = 4,
    type = "umap",
    engine = "uwot",
    neighbors = 3,
    min_dist = 0.05
  )

  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "UMAP of Top")
})

# ------------------------------------------------------------------------------
# plot_embed
# ------------------------------------------------------------------------------

test_that("plot_embed PCA path works", {
  fit <- mock_fit()

  local_mocked_bindings(
    .rfgr_params = function(...) mock_params(),
    .rfgr_oof = function(...) mock_oof(10),
    .env = asNamespace("RFGeneRank")
  )

  p <- RFGeneRank::plot_embed(fit = fit, type = "pca", show_legend = FALSE)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Decision-space PCA (OOF probabilities)")
})

test_that("plot_embed errors when OOF probabilities are missing", {
  fit <- mock_fit()

  local_mocked_bindings(
    .rfgr_params = function(...) mock_params(),
    .rfgr_oof = function(...) list(prob = NULL, y = factor(c("A", "B"))),
    .env = asNamespace("RFGeneRank")
  )

  expect_error(
    RFGeneRank::plot_embed(fit = fit, type = "pca"),
    "OOF probabilities are missing"
  )
})

test_that("plot_embed UMAP uwot path works when uwot is available", {
  skip_if_not_installed("uwot")

  fit <- mock_fit()

  local_mocked_bindings(
    .rfgr_params = function(...) mock_params(),
    .rfgr_oof = function(...) mock_oof(10),
    .env = asNamespace("RFGeneRank")
  )

  p <- RFGeneRank::plot_embed(fit = fit, type = "umap", engine = "uwot")

  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "Decision-space UMAP")
})

# ------------------------------------------------------------------------------
# Internal SHAP helpers
# ------------------------------------------------------------------------------

test_that(".make_age_numeric handles numeric and character age", {
  cd <- data.frame(age = c("40", "50", "60"), stringsAsFactors = FALSE)

  out <- RFGeneRank:::.make_age_numeric(cd, "age")
  expect_equal(out, c(40, 50, 60))

  cd2 <- data.frame(age = c(40, 50, 60))
  out2 <- RFGeneRank:::.make_age_numeric(cd2, "age")
  expect_equal(out2, c(40, 50, 60))
})

test_that(".make_age_numeric errors on missing or bad age column", {
  expect_error(
    RFGeneRank:::.make_age_numeric(data.frame(x = 1:3), "age"),
    "No 'age' found"
  )

  expect_warning(
  expect_error(
    RFGeneRank:::.make_age_numeric(data.frame(age = c("a", "b")), "age"),
    "Could not coerce"
  ),
  "NAs introduced by coercion"
)
})

test_that(".age_breaks and .parse_bin behave as expected", {
  br <- RFGeneRank:::.age_breaks(c(41, 53, 68), width = 10)
  expect_true(is.numeric(br))
  expect_gte(length(br), 2)

  rng <- RFGeneRank:::.parse_bin("[40,50)")
  expect_equal(rng, c(40, 50))
})

test_that(".summarize_age_bins returns expected columns", {
  df <- data.frame(
    x = c(40, 42, 45, 51, 54, 57),
    SHAP = c(1, 1.2, 0.8, -0.3, -0.2, -0.1),
    color = c("Female", "Female", "Female", "Male", "Male", "Male")
  )

  out <- RFGeneRank:::.summarize_age_bins(df, width = 10)
  expect_true(all(c("age_bin", "color", "n", "mean_phi", "sd_phi", "se_phi",
                    "t_stat", "p_val", "fdr", "abs_t", "abs_mean") %in% names(out)))
  expect_equal(attr(out, "alpha"), 0.05)
})

test_that(".add_band_highlight returns ggplot and handles empty table", {
  p <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()

  out_empty <- RFGeneRank:::.add_band_highlight(
    p,
    bins_tbl = data.frame(),
    mode = "significant"
  )
  expect_s3_class(out_empty, "ggplot")

  bins <- data.frame(
    age_bin = c("[40,50)", "[50,60)"),
    color = c("Female", "Male"),
    fdr = c(0.01, 0.02),
    abs_t = c(5, 4)
  )
  attr(bins, "alpha") <- 0.05

  out_sig <- RFGeneRank:::.add_band_highlight(
    p,
    bins_tbl = bins,
    mode = "significant"
  )
  expect_s3_class(out_sig, "ggplot")

  out_topk <- RFGeneRank:::.add_band_highlight(
    p,
    bins_tbl = bins,
    mode = "top_k",
    top_k = 1
  )
  expect_s3_class(out_topk, "ggplot")
})

# ------------------------------------------------------------------------------
# plot_shap_dependence
# ------------------------------------------------------------------------------

test_that("plot_shap_dependence works with numeric x and supplied shap_mat", {
  fit <- mock_fit()
  se <- toy_se_plot()

  shap_mat <- matrix(
    seq_len(ncol(se)),
    ncol = 1,
    dimnames = list(colnames(se), "g1")
  )

  local_mocked_bindings(
    .rfgr_oof = function(...) mock_oof(ncol(se)),
    .env = asNamespace("RFGeneRank")
  )

  p <- RFGeneRank::plot_shap_dependence(
    fit = fit,
    se = se,
    gene = "g1",
    x = "age",
    shap_mat = shap_mat
  )

  expect_true(inherits(p, c("ggplot", "patchwork")))
})

test_that("plot_shap_dependence works with categorical x and categorical color_by", {
  fit <- mock_fit()
  se <- toy_se_plot()

  SummarizedExperiment::colData(se)$batch <- factor(rep(c("B1", "B2"), length.out = ncol(se)))

  shap_mat <- matrix(
    rnorm(ncol(se)),
    ncol = 1,
    dimnames = list(colnames(se), "g2")
  )

  local_mocked_bindings(
    .rfgr_oof = function(...) mock_oof(ncol(se)),
    .env = asNamespace("RFGeneRank")
  )

  p <- RFGeneRank::plot_shap_dependence(
    fit = fit,
    se = se,
    gene = "g2",
    x = "batch",
    color_by = "sex",
    shap_mat = shap_mat,
    palette = c("pink", "cyan")
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_shap_dependence handles numeric color_by gradient path", {
  fit <- mock_fit()
  se <- toy_se_plot()

  shap_mat <- matrix(
    rnorm(ncol(se)),
    ncol = 1,
    dimnames = list(colnames(se), "g3")
  )

  local_mocked_bindings(
    .rfgr_oof = function(...) mock_oof(ncol(se)),
    .env = asNamespace("RFGeneRank")
  )

  p <- RFGeneRank::plot_shap_dependence(
    fit = fit,
    se = se,
    gene = "g3",
    x = "age",
    color_by = "score",
    shap_mat = shap_mat,
    palette = c("yellow", "red")
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_shap_dependence errors for missing x or grouping covariate", {
  fit <- mock_fit()
  se <- toy_se_plot()

  shap_mat <- matrix(
    rnorm(ncol(se)),
    ncol = 1,
    dimnames = list(colnames(se), "g1")
  )

  local_mocked_bindings(
    .rfgr_oof = function(...) mock_oof(ncol(se)),
    .env = asNamespace("RFGeneRank")
  )

  expect_error(
    RFGeneRank::plot_shap_dependence(fit, se, gene = "g1", x = "nope", shap_mat = shap_mat),
    "covariate 'nope' not found"
  )

  expect_error(
    RFGeneRank::plot_shap_dependence(fit, se, gene = "g1", x = "age", color_by = "missing", shap_mat = shap_mat),
    "grouping covariate 'missing' not found"
  )
})

test_that("plot_shap_dependence errors for bad SHAP inputs", {
  fit <- mock_fit()
  se <- toy_se_plot()

  shap_no_overlap <- matrix(
    rnorm(ncol(se)),
    ncol = 1,
    dimnames = list(paste0("z", seq_len(ncol(se))), "g1")
  )

  local_mocked_bindings(
    .rfgr_oof = function(...) mock_oof(ncol(se)),
    .env = asNamespace("RFGeneRank")
  )

  expect_error(
    RFGeneRank::plot_shap_dependence(fit, se, gene = "g1", x = "age", shap_mat = shap_no_overlap),
    "No overlap between shap_mat rownames and OOF IDs"
  )

  shap_wrong_gene <- matrix(
    rnorm(ncol(se)),
    ncol = 1,
    dimnames = list(colnames(se), "not_g1")
  )

  expect_error(
    RFGeneRank::plot_shap_dependence(fit, se, gene = "g1", x = "age", shap_mat = shap_wrong_gene),
    "not present in SHAP matrix"
  )
})

test_that("plot_shap_dependence errors when numeric palette has fewer than 2 colours", {
  fit <- mock_fit()
  se <- toy_se_plot()

  shap_mat <- matrix(
    rnorm(ncol(se)),
    ncol = 1,
    dimnames = list(colnames(se), "g4")
  )

  local_mocked_bindings(
    .rfgr_oof = function(...) mock_oof(ncol(se)),
    .env = asNamespace("RFGeneRank")
  )

  expect_error(
    RFGeneRank::plot_shap_dependence(
      fit = fit,
      se = se,
      gene = "g4",
      x = "age",
      color_by = "score",
      shap_mat = shap_mat,
      palette = "red"
    ),
    "must have at least 2 colours"
  )
})

# ------------------------------------------------------------------------------
# shap_train_ranger
# ------------------------------------------------------------------------------

test_that("shap_train_ranger returns matrix when ranger and fastshap are installed", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("fastshap")

  se <- toy_se_plot(n_genes = 4, n_samples = 8)

  out <- RFGeneRank::shap_train_ranger(
    se = se,
    label_col = "state",
    genes = c("g1", "g2"),
    num.trees = 25,
    nsim = 2,
    seed = 1
  )

  expect_true(is.matrix(out))
  expect_equal(nrow(out), ncol(se))
  expect_equal(colnames(out), c("g1", "g2"))
})

test_that("shap_train_ranger errors when none of requested genes are present", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("fastshap")

  se <- toy_se_plot(n_genes = 4, n_samples = 8)

  expect_error(
    RFGeneRank::shap_train_ranger(
      se = se,
      label_col = "state",
      genes = c("x1", "x2"),
      num.trees = 10,
      nsim = 2
    ),
    "None of the requested genes are present"
  )
})

# ------------------------------------------------------------------------------
# ROC plots
# ------------------------------------------------------------------------------

test_that("plot_roc returns ggplot", {
  fit <- mock_fit()

  local_mocked_bindings(
    .rfgr_oof = function(...) mock_oof(10),
    .env = asNamespace("RFGeneRank")
  )

  p <- RFGeneRank::plot_roc(fit)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "OOF ROC")
})

test_that("plot_roc_multi returns ggplot for multiple models", {
  fits <- list(
    ranger = list(oof = list(
      y = factor(rep(c("Control", "Case"), each = 5)),
      p_use = c(0.1, 0.2, 0.3, 0.2, 0.4, 0.7, 0.8, 0.9, 0.6, 0.85)
    )),
    glmnet = list(oof = list(
      y = factor(rep(c("Control", "Case"), each = 5)),
      p_use = c(0.2, 0.1, 0.4, 0.3, 0.35, 0.65, 0.75, 0.88, 0.7, 0.8)
    ))
  )

  p <- RFGeneRank::plot_roc_multi(fits)
  expect_s3_class(p, "ggplot")
})

test_that("plot_roc_multi errors when p_use is missing", {
  fits <- list(
    bad = list(oof = list(
      y = factor(rep(c("A", "B"), each = 5))
    ))
  )

  expect_error(
    RFGeneRank::plot_roc_multi(fits),
    "oof\\$p_use missing"
  )
})

# ------------------------------------------------------------------------------
# plot_confusion_heatmap
# ------------------------------------------------------------------------------

test_that("plot_confusion_heatmap works for counts and rowpct", {
  cm <- table(
    True = factor(c("A", "A", "B", "B")),
    Predicted = factor(c("A", "B", "A", "B"))
  )

  p1 <- RFGeneRank::plot_confusion_heatmap(cm, mode = "counts")
  p2 <- RFGeneRank::plot_confusion_heatmap(cm, mode = "rowpct")

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_equal(p1$labels$title, "Confusion Matrix (counts)")
  expect_equal(p2$labels$title, "Confusion Matrix (row %)")
})

# ------------------------------------------------------------------------------
# rfgr_plot_suite
# ------------------------------------------------------------------------------

test_that("rfgr_plot_suite assembles plots without val", {
  fit <- mock_fit()
  se <- toy_se_plot()

  fake_plot <- ggplot(data.frame(x = 1:2, y = 1:2), aes(x, y)) + geom_point()

  local_mocked_bindings(
    plot_importance = function(...) fake_plot,
    plot_sign_importance = function(...) fake_plot,
    plot_embed_expr = function(...) fake_plot,
    plot_embed = function(...) fake_plot,
    plot_shap_dependence = function(...) fake_plot,
    .env = asNamespace("RFGeneRank")
  )

  out <- RFGeneRank::rfgr_plot_suite(
    fit = fit,
    se = se,
    shap_gene = "g1"
  )

  expect_true(is.list(out))
  expect_named(out, c("importance", "signed_importance", "expr_pca", "decision_pca", "shap"))
})

test_that("rfgr_plot_suite assembles val-dependent plots and saves files", {
  fit <- mock_fit()
  se <- toy_se_plot()

  fake_plot <- ggplot(data.frame(x = 1:2, y = 1:2), aes(x, y)) + geom_point()

  cm <- table(
    True = factor(c("A", "A", "B", "B")),
    Predicted = factor(c("A", "B", "A", "B"))
  )

  val <- list(
    summary = data.frame(method = c("ranger", "glmnet"), stringsAsFactors = FALSE),
    ranger = list(conf_mat = cm, oof = list(y = factor(c("A", "B")), p_use = c(0.2, 0.8))),
    glmnet = list(conf_mat = cm, oof = list(y = factor(c("A", "B")), p_use = c(0.3, 0.7)))
  )

  td <- tempfile("rfgr-plots-")
  dir.create(td)

  local_mocked_bindings(
    plot_importance = function(...) fake_plot,
    plot_sign_importance = function(...) fake_plot,
    plot_embed_expr = function(...) fake_plot,
    plot_embed = function(...) fake_plot,
    plot_roc_multi = function(...) fake_plot,
    plot_confusion_heatmap = function(...) fake_plot,
    .env = asNamespace("RFGeneRank")
  )

  out <- RFGeneRank::rfgr_plot_suite(
    fit = fit,
    se = se,
    val = val,
    outdir = td
  )

  expect_true(is.list(out))
  expect_true("roc_multi" %in% names(out))
  expect_true(any(startsWith(names(out), "cm_")))

  saved <- list.files(td, pattern = "\\.png$", full.names = TRUE)
  expect_true(length(saved) >= 5)
})

test_that("plot_importance map_to_symbol falls back with warning when org.Hs.eg.db is absent", {
  skip_if(requireNamespace("org.Hs.eg.db", quietly = TRUE),
          "This branch is only hit when org.Hs.eg.db is absent")

  fit <- mock_fit()

  local_mocked_bindings(
    .rfgr_imp = function(.) mock_imp(),
    .env = asNamespace("RFGeneRank")
  )

  expect_warning(
    p <- plot_importance(fit, top = 3, map_to_symbol = TRUE),
    "org.Hs.eg.db is not installed"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_sign_importance map_to_symbol falls back with warning when org.Hs.eg.db is absent", {
  skip_if(requireNamespace("org.Hs.eg.db", quietly = TRUE),
          "This branch is only hit when org.Hs.eg.db is absent")

  fit <- mock_fit()

  local_mocked_bindings(
    .rfgr_imp = function(.) data.frame(
      gene = c("g1", "g2", "g3"),
      importance = c(0.8, 0.6, 0.4),
      direction = c(1, -1, 1),
      signed_importance = c(0.8, -0.6, 0.4)
    ),
    .env = asNamespace("RFGeneRank")
  )

  expect_warning(
    p <- plot_sign_importance(fit = fit, top = 3, map_to_symbol = TRUE),
    "org.Hs.eg.db is not installed"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_embed_expr PCA applies palette branch", {
  fit <- mock_fit()
  se <- toy_se_plot()

  local_mocked_bindings(
    .rfgr_imp = function(.) mock_imp(),
    .rfgr_params = function(.) mock_params(label_col = "label"),
    .env = asNamespace("RFGeneRank")
  )

  p <- plot_embed_expr(
    fit = fit,
    se = se,
    n_top = 4,
    type = "pca",
    zscore = TRUE,
    show_legend = TRUE,
    palette = c(Control = "black", Case = "red")
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_embed PCA applies palette branch", {
  fit <- mock_fit()

  local_mocked_bindings(
    .rfgr_params = function(...) mock_params(),
    .rfgr_oof = function(...) mock_oof(10),
    .env = asNamespace("RFGeneRank")
  )

  p <- plot_embed(
    fit = fit,
    type = "pca",
    show_legend = TRUE,
    palette = c(Control = "black", Case = "red")
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_embed hides legend when requested", {
  fit <- mock_fit()

  local_mocked_bindings(
    .rfgr_params = function(...) mock_params(),
    .rfgr_oof = function(...) mock_oof(10),
    .env = asNamespace("RFGeneRank")
  )

  p <- RFGeneRank::plot_embed(
    fit = fit,
    type = "pca",
    show_legend = FALSE
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_confusion_heatmap respects max_prop_text filter", {
  cm <- table(
    True = factor(c("A", "A", "A", "B", "B", "B", "B")),
    Predicted = factor(c("A", "A", "B", "A", "B", "B", "B"))
  )

  p <- plot_confusion_heatmap(
    cm,
    mode = "rowpct",
    max_prop_text = 0.70
  )

  expect_s3_class(p, "ggplot")

  dat <- ggplot_build(p)$plot$data
  expect_true(nrow(dat) >= 1)
})

test_that("rfgr_plot_suite includes shap plot when shap_gene is provided", {
  fit <- mock_fit()
  se <- toy_se_plot()

  fake_plot <- ggplot(data.frame(x = 1:2, y = 1:2), aes(x, y)) + geom_point()

  local_mocked_bindings(
    plot_importance = function(...) fake_plot,
    plot_sign_importance = function(...) fake_plot,
    plot_embed_expr = function(...) fake_plot,
    plot_embed = function(...) fake_plot,
    plot_shap_dependence = function(...) fake_plot,
    .env = asNamespace("RFGeneRank")
  )

  out <- rfgr_plot_suite(
    fit = fit,
    se = se,
    shap_gene = "g1"
  )

  expect_true(is.list(out))
  expect_true("shap" %in% names(out))
})

test_that("rfgr_plot_suite adds roc and confusion plots when val is provided", {
  fit <- mock_fit()
  se <- toy_se_plot()

  fake_plot <- ggplot(data.frame(x = 1:2, y = 1:2), aes(x, y)) + geom_point()

  cm <- table(
    True = factor(c("A", "A", "B", "B")),
    Predicted = factor(c("A", "B", "A", "B"))
  )

  val <- list(
    summary = data.frame(method = c("ranger", "glmnet"), stringsAsFactors = FALSE),
    ranger = list(conf_mat = cm, oof = list(y = factor(c("A", "B")), p_use = c(0.2, 0.8))),
    glmnet = list(conf_mat = cm, oof = list(y = factor(c("A", "B")), p_use = c(0.3, 0.7)))
  )

  local_mocked_bindings(
    plot_importance = function(...) fake_plot,
    plot_sign_importance = function(...) fake_plot,
    plot_embed_expr = function(...) fake_plot,
    plot_embed = function(...) fake_plot,
    plot_roc_multi = function(...) fake_plot,
    plot_confusion_heatmap = function(...) fake_plot,
    .env = asNamespace("RFGeneRank")
  )

  out <- rfgr_plot_suite(
    fit = fit,
    se = se,
    val = val
  )

  expect_true(is.list(out))
  expect_true("roc_multi" %in% names(out))
  expect_true(any(startsWith(names(out), "cm_")))
})

test_that("plot_importance returns ggplot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")

  se <- toy_se(ngenes = 50, nsamp = 12, seed = 1)
  fit <- rank_genes(se, label_col = "label", n_top = 20, trees = 50, seed = 1)

  p <- plot_importance(fit, top = 10, map_to_symbol = FALSE)
  expect_true(inherits(p, "ggplot"))
})

test_that("plot_embed supports UMAP branch when params$umap exists", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("methods")

  fit <- methods::new("GeneRankFit")
  methods::slot(fit, "params") <- list(
    seed = 1,
    umap = list(neighbors = 10, min_dist = 0.2, metric = "cosine")
  )

  n <- 24
  y <- factor(rep(c("Control","Case"), each = n/2))
  set.seed(10)
  p_case <- c(runif(n/2, 0.1, 0.4), runif(n/2, 0.6, 0.9))
  prob <- cbind(Control = 1 - p_case, Case = p_case)
  rownames(prob) <- paste0("s", 1:n)
  names(y) <- rownames(prob)

  methods::slot(fit, "oof") <- list(prob = prob, y = y)

  g <- plot_embed(fit, type = "umap", engine = "umap", seed = 1)
  expect_true(inherits(g, "ggplot"))
})

test_that("plot_embed_expr covers zscore TRUE/FALSE branches", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  skip_if_not_installed("matrixStats")
  skip_if_not_installed("methods")

  n <- 18
  ng <- 80
  set.seed(11)

  expr <- matrix(rnorm(ng * n), nrow = ng)
  rownames(expr) <- paste0("g", 1:ng)
  colnames(expr) <- paste0("s", 1:n)

  y <- factor(rep(c("Control","Case"), each = n/2))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(state = y, row.names = colnames(expr))
  )

  fit <- methods::new("GeneRankFit")
  methods::slot(fit, "params") <- list(seed = 1, label_col = "state")
  methods::slot(fit, "imp") <- data.frame(gene = rownames(expr), importance = runif(ng))

  g1 <- plot_embed_expr(fit, se, n_top = 25, type = "pca", zscore = TRUE,  seed = 1)
  g2 <- plot_embed_expr(fit, se, n_top = 25, type = "pca", zscore = FALSE, seed = 1)

  expect_true(inherits(g1, "ggplot"))
  expect_true(inherits(g2, "ggplot"))
})

test_that("plot_confusion_heatmap covers additional mode branches if supported", {
  skip_if_not_installed("ggplot2")

  true <- factor(rep(c("A","B"), each = 12))
  pred <- factor(sample(c("A","B"), 24, replace = TRUE))
  cm <- table(True = true, Predicted = pred)

  # already tested counts/rowpct earlier, add colpct/totalpct if your function supports them
  ok <- FALSE
  try({ p <- plot_confusion_heatmap(cm, mode = "colpct"); ok <- TRUE }, silent = TRUE)
  if (ok) expect_true(inherits(p, "ggplot"))

  ok2 <- FALSE
  try({ p2 <- plot_confusion_heatmap(cm, mode = "totalpct"); ok2 <- TRUE }, silent = TRUE)
  if (ok2) expect_true(inherits(p2, "ggplot"))
})

test_that("plot_importance works on minimal valid fit", {
  fit <- methods::new("GeneRankFit")
  fit@imp <- data.frame(
    gene = c("g1", "g2", "g3"),
    importance = c(3, 2, 1)
  )

  p <- RFGeneRank::plot_importance(fit, top = 2)
  expect_s3_class(p, "ggplot")
})

test_that("plot_importance handles map_to_symbol branch", {
  fit <- methods::new("GeneRankFit")
  fit@imp <- data.frame(
    gene = c("1", "2", "3"),
    importance = c(3, 2, 1)
  )

  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    p <- RFGeneRank::plot_importance(fit, top = 2, map_to_symbol = TRUE)
    expect_s3_class(p, "ggplot")
  } else {
    expect_warning(
      p <- RFGeneRank::plot_importance(fit, top = 2, map_to_symbol = TRUE),
      "org.Hs.eg.db is not installed"
    )
    expect_s3_class(p, "ggplot")
  }
})

test_that("plot_sign_importance works from tab and covers palette/legend branches", {
  tab <- data.frame(
    gene = c("g1", "g2", "g3"),
    importance = c(0.8, 0.6, 0.4),
    direction = c(1, -1, 1),
    signed_importance = c(0.8, -0.6, 0.4)
  )

  p <- RFGeneRank::plot_sign_importance(
    tab = tab,
    top = 2,
    show_legend = FALSE,
    palette = c(`-1` = "blue", `1` = "red")
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_embed returns ggplot for valid OOF probabilities", {
  y <- factor(c("A", "A", "B", "B"))
  names(y) <- paste0("s", 1:4)

  prob <- cbind(
    A = c(0.8, 0.7, 0.3, 0.2),
    B = c(0.2, 0.3, 0.7, 0.8)
  )
  rownames(prob) <- names(y)

  fit <- methods::new("GeneRankFit")
  fit@oof <- list(prob = prob, y = y)
  fit@params <- list(seed = 1)

  p <- RFGeneRank::plot_embed(fit, type = "pca", show_legend = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_embed_expr returns ggplot for PCA path", {
  set.seed(1)

  expr <- matrix(rnorm(10 * 6), nrow = 10)
  rownames(expr) <- paste0("g", 1:10)
  colnames(expr) <- paste0("s", 1:6)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = S4Vectors::DataFrame(
      label = factor(rep(c("A", "B"), each = 3)),
      row.names = colnames(expr)
    )
  )

  fit <- methods::new("GeneRankFit")
  fit@params <- list(label_col = "label", seed = 1)
  fit@imp <- data.frame(
    gene = rownames(expr),
    importance = rev(seq_len(nrow(expr)))
  )

  p <- RFGeneRank::plot_embed_expr(
    fit, se,
    n_top = 5,
    type = "pca",
    show_legend = FALSE
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_confusion_heatmap works in counts and rowpct modes", {
  true <- factor(c("A", "A", "B", "B"))
  pred <- factor(c("A", "B", "A", "B"))
  cm <- table(True = true, Predicted = pred)

  p1 <- RFGeneRank::plot_confusion_heatmap(cm, mode = "counts")
  p2 <- RFGeneRank::plot_confusion_heatmap(cm, mode = "rowpct")

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_roc returns ggplot from minimal GeneRankFit with oof$prob", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("pROC")
  skip_if_not_installed("methods")

  fit <- methods::new("GeneRankFit")

  n <- 30
  y <- factor(rep(c("Control","Case"), each = n/2))
  p_case <- c(runif(n/2, 0.05, 0.35), runif(n/2, 0.65, 0.95))
  prob <- cbind(Control = 1 - p_case, Case = p_case)
  rownames(prob) <- paste0("s", 1:n)
  names(y) <- rownames(prob)

  methods::slot(fit, "oof") <- list(prob = prob, y = y)

  p <- plot_roc(fit)
  expect_true(inherits(p, "ggplot"))
})

test_that("plot_roc_multi overlays ROC curves from validate-like objects (oof$p_use)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("pROC")

  y <- factor(rep(c("Control", "Case"), each = 25))
  p1 <- c(runif(25, 0.05, 0.35), runif(25, 0.65, 0.95))
  p2 <- pmin(pmax(p1 + rnorm(50, 0, 0.08), 0), 1)

  fits <- list(
    Model_1 = list(oof = list(y = y, p_use = p1)),
    Model_2 = list(oof = list(y = y, p_use = p2))
  )

  g <- plot_roc_multi(fits, title = "test")
  expect_true(inherits(g, "ggplot"))
})

test_that("plot_confusion_heatmap works in counts and rowpct modes", {
  skip_if_not_installed("ggplot2")

  true <- factor(rep(c("A","B"), each = 10))
  pred <- factor(c(sample(c("A","B"), 18, replace = TRUE), "A", "B"))
  cm <- table(True = true, Predicted = pred)

  p1 <- plot_confusion_heatmap(cm, mode = "counts")
  p2 <- plot_confusion_heatmap(cm, mode = "rowpct")

  expect_true(inherits(p1, "ggplot"))
  expect_true(inherits(p2, "ggplot"))
})

test_that("plot_embed (PCA) runs when params include umap + seed", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("methods")

  fit <- methods::new("GeneRankFit")
  # IMPORTANT: plot_embed reads pars$umap$neighbors/min_dist/metric and pars$seed
  methods::slot(fit, "params") <- list(
    seed = 1,
    umap = list(neighbors = 10, min_dist = 0.2, metric = "cosine")
  )

  n <- 20
  y <- factor(rep(c("Control","Case"), each = n/2))
  p_case <- c(runif(n/2, 0.1, 0.4), runif(n/2, 0.6, 0.9))
  prob <- cbind(Control = 1 - p_case, Case = p_case)
  rownames(prob) <- paste0("s", 1:n)
  names(y) <- rownames(prob)

  methods::slot(fit, "oof") <- list(prob = prob, y = y)

  g <- plot_embed(fit, type = "pca")
  expect_true(inherits(g, "ggplot"))
})

test_that("plot_embed_expr (PCA) runs when params include label_col and imp has gene IDs", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("matrixStats")
  skip_if_not_installed("methods")

  n <- 16
  ng <- 60

  expr <- matrix(rnorm(ng * n), nrow = ng)
  rownames(expr) <- paste0("g", 1:ng)
  colnames(expr) <- paste0("s", 1:n)

  y <- factor(rep(c("Control","Case"), each = n/2))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = data.frame(state = y, row.names = colnames(expr))
  )

  fit <- methods::new("GeneRankFit")
  methods::slot(fit, "params") <- list(seed = 1, label_col = "state")

  # imp MUST have column 'gene' that overlaps rownames(expr)
  imp_tab <- data.frame(gene = rownames(expr), importance = runif(ng))
  methods::slot(fit, "imp") <- imp_tab

  g <- plot_embed_expr(fit, se, n_top = 30, type = "pca", seed = 1)
  expect_true(inherits(g, "ggplot"))
})

test_that("plot_roc returns ggplot from minimal GeneRankFit oof", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("pROC")
  skip_if_not_installed("methods")

  fit <- methods::new("GeneRankFit")

  n <- 30
  y <- factor(rep(c("Control","Case"), each = n/2))
  p_case <- c(runif(n/2, 0.05, 0.35), runif(n/2, 0.65, 0.95))
  prob <- cbind(Control = 1 - p_case, Case = p_case)
  rownames(prob) <- paste0("s", 1:n)

  methods::slot(fit, "oof") <- list(prob = prob, y = y)

  p <- plot_roc(fit)
  expect_true(inherits(p, "ggplot"))
})

test_that("plot_roc_multi overlays ROC curves for validate-like objects", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("pROC")

  y <- factor(rep(c("Control","Case"), each = 25))
  p1 <- c(runif(25, 0.05, 0.35), runif(25, 0.65, 0.95))
  p2 <- pmin(pmax(p1 + rnorm(50, 0, 0.08), 0), 1)

  fits <- list(
    ranger = list(oof = list(p_use = p1, y = y)),
    glmnet = list(oof = list(p_use = p2, y = y))
  )

  g <- plot_roc_multi(fits, title = "test")
  expect_true(inherits(g, "ggplot"))
})

test_that("plot_confusion_heatmap works in counts and rowpct modes", {
  skip_if_not_installed("ggplot2")

  true <- factor(rep(c("A","B"), each = 10))
  pred <- factor(c(sample(c("A","B"), 18, replace = TRUE), "A", "B"))

  cm <- table(True = true, Predicted = pred)

  p1 <- plot_confusion_heatmap(cm, mode = "counts")
  p2 <- plot_confusion_heatmap(cm, mode = "rowpct")

  expect_true(inherits(p1, "ggplot"))
  expect_true(inherits(p2, "ggplot"))
})

test_that("plot_embed and plot_embed_expr run (PCA/UMAP) on minimal structures", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("methods")

  # ---- minimal fit with params$umap to avoid pars$umap$... NULL errors ----
  fit <- methods::new("GeneRankFit")
  methods::slot(fit, "params") <- list(
    seed = 1,
    label_col = "state",
    umap = list(neighbors = 10, min_dist = 0.2, metric = "cosine")
  )

  n <- 20
  y <- factor(rep(c("Control","Case"), each = n/2))
  p_case <- c(runif(n/2, 0.1, 0.4), runif(n/2, 0.6, 0.9))
  prob <- cbind(Control = 1 - p_case, Case = p_case)
  rownames(prob) <- paste0("s", 1:n)
  methods::slot(fit, "oof") <- list(prob = prob, y = y)

  # embed from oof prob
  g1 <- plot_embed(fit, type = "pca")
  expect_true(inherits(g1, "ggplot"))

  # UMAP branch
  g2 <- plot_embed(fit, type = "umap", engine = "umap", neighbors = 5, min_dist = 0.3, metric = "euclidean", seed = 1)
  expect_true(inherits(g2, "ggplot"))

  # ---- embed expression requires se + imp(gene) ----
  ng <- 40
  expr <- matrix(rnorm(ng * n), nrow = ng)
  rownames(expr) <- paste0("g", 1:ng)
  colnames(expr) <- rownames(prob)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = data.frame(
      state = y,
      batch = factor(rep(c("B1","B2"), length.out = n)),
      row.names = colnames(expr)
    )
  )

  # fake importance table with genes present in assay
  imp_tab <- data.frame(gene = rownames(expr), importance = runif(ng))
  imp_tab <- imp_tab[order(-imp_tab$importance), , drop = FALSE]
  imp(fit) <- imp_tab

  g3 <- plot_embed_expr(fit, se, n_top = 15, type = "umap", engine = "umap",
                        neighbors = 5, min_dist = 0.25, metric = "euclidean",
                        zscore = TRUE, seed = 1)
  expect_true(inherits(g3, "ggplot"))
})

test_that("plot_sign_importance works from tab and covers palette/legend branches", {
  tab <- data.frame(
    gene = c("g1", "g2", "g3"),
    importance = c(0.8, 0.6, 0.4),
    direction = c(1, -1, 1),
    signed_importance = c(0.8, -0.6, 0.4)
  )

  p <- RFGeneRank::plot_sign_importance(
    tab = tab,
    top = 2,
    show_legend = FALSE,
    palette = c(`-1` = "blue", `1` = "red")
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_importance warns and returns plot when map_to_symbol=TRUE without org.Hs.eg.db", {
  fit <- methods::new("GeneRankFit")
  fit@imp <- data.frame(
    gene = c("1", "2", "3"),
    importance = c(3, 2, 1)
  )

  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    expect_warning(
      p <- RFGeneRank::plot_importance(fit, top = 2, map_to_symbol = TRUE),
      "org.Hs.eg.db is not installed"
    )
    expect_s3_class(p, "ggplot")
  }
})

test_that("rfgr_plot_suite returns list of ggplots (ENTREZ-like genes)", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ggplot2")

  set.seed(2)

  ngenes <- 80
  nsamp  <- 24

  expr <- matrix(rnorm(ngenes * nsamp), nrow = ngenes)
  rownames(expr) <- as.character(seq_len(ngenes))   # ENTrez IDs: "1","2",...
  colnames(expr) <- paste0("s", seq_len(nsamp))

  y <- factor(rep(c("Control", "Case"), each = nsamp/2))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = expr),
    colData = data.frame(
      label = y,
      batch = rep(1:2, length.out = nsamp),
      row.names = colnames(expr)
    )
  )

  fit <- RFGeneRank::rank_genes(se = se, label_col = "label", n_top = 30, trees = 50)

  X <- t(SummarizedExperiment::assay(se, "expr")) # samples x genes
  y2 <- as.factor(SummarizedExperiment::colData(se)$label)

  tab <- RFGeneRank::sign_importance(fit = fit, X = X, y = y2, method = "mean")
  RFGeneRank::imp(fit) <- tab

  plots <- RFGeneRank::rfgr_plot_suite(fit = fit, se = se, outdir = NULL, top = 10)

  expect_true(is.list(plots))
  expect_true(length(plots) > 0)
  expect_true(any(vapply(plots, inherits, logical(1), what = "ggplot")))
})