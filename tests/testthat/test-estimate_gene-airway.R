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
    draws_warmup = 150,
    draws_sampling = 150,
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

  p <- suppressMessages(plot_boxplot(fit, factor = "dex", number_of_draws = 20))
  expect_s3_class(p, "ggplot")
  p_adj <- plot_boxplot(
    fit,
    factor = "dex",
    remove_unwanted_effects = TRUE,
    number_of_draws = 20
  )
  expect_s3_class(p_adj, "ggplot")
})

test_that("ZINB with a dispersion-derived shape prior fits one airway gene", {
  skip_on_cran()
  skip_if_no_cmdstan()
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")

  se <- airway_se(n_genes = 150)
  se <- estimate_dispersion_prior(
    se,
    formula_abundance = ~ dex,
    method = "degrees_freedom"
  )
  se$offset <- log(colSums(SummarizedExperiment::assay(se, "counts")))
  se_gene <- se["ENSG00000120129", , drop = FALSE]

  fit <- estimate_gene(
    se_gene,
    formula_abundance = ~ dex,
    offset = "offset",
    dispersion_prior_log_mean = "dispersion_prior_log_mean",
    dispersion_prior_log_sd = "dispersion_prior_log_sd",
    chains = 1,
    draws_warmup = 100,
    draws_sampling = 150,
    cores = 1,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  expect_s3_class(fit, "brmsfit")
  expect_true("dispersion_prior_log_mean" %in% names(fit$data))
  expect_true("b_shape_Intercept" %in% brms::variables(fit))

  shape_prior <- brms::prior_summary(fit)
  shape_int <- shape_prior$prior[
    shape_prior$class == "Intercept" & shape_prior$dpar == "shape"
  ]
  expect_identical(shape_int, "student_t(3, 0, brmde_shape_scale)")
  expect_true(is.finite(brms::standata(fit)$brmde_shape_scale))
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
    draws_warmup = 100,
    draws_sampling = 150,
    cores = 1,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  expect_s3_class(fit, "brmsfit")
  expect_equal(fit$family$family, "zero_inflated_negbinomial")
  expect_true("b_dextrt" %in% brms::variables(fit))
})

test_that("hypothesis_gene reports the convergence of each contrast", {
  skip_on_cran()
  skip_if_no_cmdstan()

  se <- airway_one_gene("ENSG00000120129")

  fit <- estimate_gene(
    se,
    formula_abundance = ~ dex + (1 | cell),
    family = brms::negbinomial(),
    offset = "offset",
    chains = 2,
    draws_warmup = 150,
    draws_sampling = 150,
    cores = 1,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  hyp <- hypothesis_gene(fit, c("dextrt = 0", "Intercept"))
  expect_equal(nrow(hyp), 2L)
  expect_true(all(c("rhat", "ess_bulk", "mcse") %in% names(hyp)))
  expect_true(all(is.finite(hyp$rhat)))
  expect_true(all(hyp$ess_bulk > 0))

  # The diagnostics belong to the contrast, so the first row's are those of
  # b_dextrt itself: that hypothesis is the parameter, unchanged. Which also
  # says the chains were recovered from brms' flattened draws correctly, since
  # every one of these statistics depends on which chain a draw came from.
  draws <- posterior::subset_draws(
    posterior::as_draws_array(fit),
    variable = "b_dextrt"
  )
  reference <- posterior::summarise_draws(
    draws,
    "rhat",
    "ess_bulk",
    "mcse_median",
    "median"
  )
  expect_equal(hyp$rhat[[1]], reference$rhat)
  expect_equal(hyp$ess_bulk[[1]], reference$ess_bulk)
  # robust = TRUE by default, so `mcse` is that of the median it reports.
  expect_equal(hyp$mcse[[1]], reference$mcse_median)
  expect_equal(hyp$estimate[[1]], reference$median)

  # The two contrasts are different quantities and get their own numbers.
  expect_false(hyp$rhat[[1]] == hyp$rhat[[2]])
})
