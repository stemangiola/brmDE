test_that("estimate_gene, hypothesis_gene, and adjust_gene run on one airway gene", {
  skip_on_cran()
  skip_if_no_cmdstan()

  se <- airway_one_gene("ENSG00000120129")
  expect_equal(nrow(se), 1L)
  expect_equal(ncol(se), 8L)
  expect_equal(rownames(se), "ENSG00000120129")

  fit <- estimate_gene(
    se,
    formula_abundance = ~ dex + (1 | cell),
    family = brms::negbinomial(),
    abundance = "counts",
    offset = "offset",
    chains = 1,
    iter = 300,
    warmup = 150,
    cores = 1,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  expect_s3_class(fit, "brmsfit")
  expect_equal(nrow(fit$data), 8L)
  expect_true("dex" %in% names(fit$data))
  expect_true("offset" %in% names(fit$data))
  expect_match(
    paste(deparse(fit$formula$formula), collapse = " "),
    "offset\\s*\\(\\s*offset\\s*\\)"
  )
  expect_true("b_dextrt" %in% brms::variables(fit))

  hyp <- hypothesis_gene(fit, "dextrt = 0")
  expect_s3_class(hyp, "tbl_df")
  expect_equal(nrow(hyp), 1L)
  expect_equal(hyp$component, "fixed")
  expect_true(is.finite(hyp$estimate))

  adj <- adjust_gene(
    fit,
    nullify = "dex",
    re_formula = ~(1 | cell),
    sample_id = colnames(se)
  )
  expect_equal(nrow(adj), 8L)
  expect_equal(adj$sample_id, colnames(se))
  expect_true(all(
    c("adjusted___Estimate", "residuals___Estimate", "fitted___Estimate") %in%
      names(adj)
  ))
})

test_that("ZINB with a dispersion-derived shape prior fits one airway gene", {
  skip_on_cran()
  skip_if_no_cmdstan()
  skip_if_not_installed("edgeR")

  se <- airway_se(n_genes = 150)
  se <- estimate_dispersion(se, ~ dex, abundance = "counts")
  se$offset <- log(colSums(SummarizedExperiment::assay(se, "counts")))
  se_gene <- se["ENSG00000120129", , drop = FALSE]

  fit <- estimate_gene(
    se_gene,
    formula_abundance = ~ dex,
    offset = "offset",
    dispersion = "dispersion",
    shape_prior_df = 3,
    chains = 1,
    iter = 250,
    warmup = 100,
    cores = 1,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  expect_s3_class(fit, "brmsfit")
  expect_true("dispersion" %in% names(fit$data))
  expect_true("b_shape_Intercept" %in% brms::variables(fit))

  shape_prior <- brms::prior_summary(fit)
  shape_int <- shape_prior$prior[
    shape_prior$class == "Intercept" & shape_prior$dpar == "shape"
  ]
  expect_match(shape_int, "^student_t\\(3, 0, 0\\.")
})

test_that("ZINB with a conjugate gamma shape prior fits one airway gene", {
  skip_on_cran()
  skip_if_no_cmdstan()
  skip_if_not_installed("edgeR")

  se <- airway_se(n_genes = 150)
  se <- estimate_dispersion(se, ~ dex, abundance = "counts")
  se$offset <- log(colSums(SummarizedExperiment::assay(se, "counts")))
  se_gene <- se["ENSG00000120129", , drop = FALSE]

  fit <- estimate_gene(
    se_gene,
    formula_abundance = ~ dex,
    offset = "offset",
    dispersion = "dispersion",
    shape_prior = "gamma",
    chains = 1,
    iter = 250,
    warmup = 100,
    cores = 1,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  expect_s3_class(fit, "brmsfit")
  # shape is sampled directly, with no submodel and no shape intercept.
  expect_true("shape" %in% brms::variables(fit))
  expect_false("b_shape_Intercept" %in% brms::variables(fit))

  prior_tbl <- brms::prior_summary(fit)
  expect_match(prior_tbl$prior[prior_tbl$class == "shape"], "^gamma\\(")

  # The prior mean is edgeR's own point estimate, 1 / phi.
  phi <- SummarizedExperiment::rowData(se_gene)$dispersion
  d_eff <- SummarizedExperiment::rowData(se_gene)$dispersion_degrees_freedom
  expect_equal((d_eff / 2) / (d_eff * phi / 2), 1 / phi)
})

test_that("estimate_gene default ZINB priors work on one airway gene", {
  skip_on_cran()
  skip_if_no_cmdstan()

  se <- airway_one_gene("ENSG00000120129")

  fit <- estimate_gene(
    se,
    formula_abundance = ~ dex,
    offset = "offset",
    chains = 1,
    iter = 250,
    warmup = 100,
    cores = 1,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  expect_s3_class(fit, "brmsfit")
  expect_equal(fit$family$family, "zero_inflated_negbinomial")
  expect_true("b_dextrt" %in% brms::variables(fit))
})
