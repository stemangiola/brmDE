# A point hypothesis is a Savage-Dickey density ratio, so it needs prior draws
# and a prior that integrates to one. These tests pin down both halves: that
# the default prior set is proper for every negative binomial family, and that
# hypothesis_gene() says so rather than returning a silent NA when it is not.

coefficient_prior_string <- function(priors) {
  priors$prior[priors$class == "b" & priors$dpar == ""]
}

test_that("the default coefficient prior is one log2 fold change wide", {
  dat <- airway_one_gene_tbl()
  priors <- brmDE:::location_priors(dat, "counts", "offset")$prior
  expect_equal(coefficient_prior_string(priors), "student_t(3, 0, 0.7)")
})

test_that("coefficient_prior_scale and coefficient_prior_df reach the prior", {
  dat <- airway_one_gene_tbl()
  priors <- brmDE:::location_priors(
    dat,
    "counts",
    "offset",
    coefficient_prior_scale = 2.5,
    coefficient_prior_df = 4
  )$prior
  expect_equal(coefficient_prior_string(priors), "student_t(4, 0, 2.5)")
})

test_that("default priors cover negbinomial as well as its zero-inflated form", {
  expect_true(brmDE:::is_negbinomial_family(brms::negbinomial()))
  expect_true(brmDE:::is_negbinomial_family(brms::zero_inflated_negbinomial()))
  expect_false(brmDE:::is_negbinomial_family(stats::poisson()))

  # Without this, a negbinomial() fit keeps brms' flat prior on `b` and no
  # amount of prior sampling can produce a Savage-Dickey ratio.
  dat <- airway_one_gene_tbl()
  formula <- brmDE:::prepare_formula(
    ~ dex,
    abundance = "counts",
    offset = "offset"
  )
  priors <- brmDE:::default_gene_priors(
    data = dat,
    formula = formula,
    abundance = "counts",
    offset = "offset",
    dispersion = NULL,
    dispersion_degrees_freedom = NULL,
    shape_prior_df = 3,
    shape_prior = "student_t"
  )$prior
  expect_true("b" %in% priors$class)
  expect_true(all(nzchar(coefficient_prior_string(priors))))
})

test_that("a point hypothesis reports pH0 when priors were sampled", {
  skip_cmdstan()

  se <- airway_one_gene("ENSG00000120129")
  fit <- estimate_gene(
    se,
    formula_abundance = ~ dex,
    family = brms::negbinomial(),
    abundance = "counts",
    offset = "offset",
    sample_prior = "yes",
    chains = 2,
    iter = 1000,
    warmup = 500,
    cores = 2,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  expect_true("prior_b" %in% brms::variables(fit))

  point <- hypothesis_gene(fit, "dextrt = 0")
  expect_equal(point$test, "point")
  expect_true(is.finite(point$pH0))
  expect_true(is.finite(point$evid_ratio))

  # brms reports the Bayes factor in favour of the null; evid_ratio is the
  # odds against it, so the two are reciprocal.
  brms_bf <- brms::hypothesis(fit, "dextrt = 0")$hypothesis$Evid.Ratio
  expect_equal(point$evid_ratio, 1 / brms_bf, tolerance = 1e-6)
  expect_equal(point$pH0, brms::hypothesis(fit, "dextrt = 0")$hypothesis$Post.Prob)
})

test_that("the default fit carries no prior draws, so a point hypothesis is NA", {
  skip_cmdstan()

  se <- airway_one_gene("ENSG00000120129")
  fit <- estimate_gene(
    se,
    formula_abundance = ~ dex,
    family = brms::negbinomial(),
    abundance = "counts",
    offset = "offset",
    chains = 2,
    iter = 1000,
    warmup = 500,
    cores = 2,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  expect_false(any(startsWith(brms::variables(fit), "prior_")))
  point <- hypothesis_gene(fit, "dextrt = 0")
  expect_true(is.na(point$pH0))
  expect_true(is.na(point$evid_ratio))

  # The effect size comes from the posterior, so it survives even when the
  # probability cannot be computed.
  expect_true(is.finite(point$log2_fold_change))

  # The directional test needs only the posterior, which is why not storing
  # prior draws is the default.
  expect_true(is.finite(hypothesis_gene(fit, "dextrt")$pH0))
})

test_that("each cell intercept can be tested as its own hypothesis", {
  skip_cmdstan()

  se <- airway_one_gene("ENSG00000120129")
  fit <- estimate_gene(
    se,
    formula_abundance = ~ dex + (1 | cell),
    family = brms::negbinomial(),
    abundance = "counts",
    offset = "offset",
    chains = 1,
    iter = 400,
    warmup = 200,
    cores = 1,
    backend = "cmdstanr",
    refresh = 0,
    silent = 2
  )

  cells <- levels(fit$data$cell)
  hyp <- paste0("`", "cell[", cells, ",Intercept]`")
  out <- hypothesis_gene(fit, hyp, class = "r")

  expect_equal(nrow(out), length(cells))
  expect_true(all(is.finite(out$estimate)))
  expect_true(all(is.finite(out$pH0)))
  expect_equal(out$hypothesis, hyp)
})
