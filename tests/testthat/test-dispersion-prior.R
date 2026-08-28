test_that("check_and_install_packages is a no-op when the package is present", {
  expect_silent(brmDE:::check_and_install_packages("stats"))
})

test_that("estimate_dispersion_prior is an S4 method for SummarizedExperiment", {
  expect_error(
    estimate_dispersion_prior(data.frame(x = 1), ~ x),
    "unable to find an inherited method"
  )
})

test_that("estimate_dispersion_prior writes prior mean and SD columns", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  skip_if_not_installed("airway")

  se <- airway_se(n_genes = 40)
  se <- estimate_dispersion_prior(se, formula_abundance = ~ dex + cell)
  rd <- SummarizedExperiment::rowData(se)
  sd <- rd$dispersion_prior_log_sd
  loc <- rd$dispersion_prior_log_mean
  expect_true(all(is.finite(sd)))
  expect_true(all(sd > 0))
  expect_true(all(sd < 5))
  expect_true(all(is.finite(loc)))
  expect_true(all(loc > 0))
  expect_null(rd$dispersion_df_prior_sd)
  expect_null(rd$dispersion_effective_degrees_freedom)

  design <- stats::model.matrix(
    ~ dex + cell,
    data = as.data.frame(SummarizedExperiment::colData(se))
  )
  shared <- brmDE:::edger_shared_loglik(
    SummarizedExperiment::assay(se, "counts"),
    design
  )
  expect_equal(loc, shared$trended.dispersion)
  residual_df <- ncol(se) - ncol(design)
  ok <- is.finite(shared$prior.n)
  expect_equal(
    shared$effective_degrees_freedom[ok],
    residual_df * (1 + shared$prior.n[ok])
  )
  # Moderated: residual df of ~ dex + cell on 8 samples is 3; d0 = prior.df > 0.
  expect_true(all(shared$effective_degrees_freedom[ok] > 3))
})

test_that("estimate_dispersion_prior method degrees_freedom writes the trigamma SD", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  skip_if_not_installed("airway")

  se <- airway_se(n_genes = 40)
  se <- estimate_dispersion_prior(
    se,
    formula_abundance = ~ dex + cell,
    method = "degrees_freedom"
  )
  rd <- SummarizedExperiment::rowData(se)
  design <- stats::model.matrix(
    ~ dex + cell,
    data = as.data.frame(SummarizedExperiment::colData(se))
  )
  shared <- brmDE:::edger_shared_loglik(
    SummarizedExperiment::assay(se, "counts"),
    design
  )
  expect_equal(
    rd$dispersion_prior_log_sd,
    brmDE:::dispersion_log_sd_from_degrees_freedom(shared$effective_degrees_freedom)
  )
  expect_equal(rd$dispersion_prior_log_mean, shared$trended.dispersion)
})

test_that("estimate_dispersion_prior rejects duplicate column names", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  skip_if_not_installed("airway")

  se <- airway_se(n_genes = 8)
  expect_error(
    estimate_dispersion_prior(
      se,
      formula_abundance = ~ dex + cell,
      log_sd_column = "x",
      mean_column = "x"
    ),
    "must name different columns"
  )
})

test_that("the dispersion grid defaults to spline.pts in [-20, 10] with 61 points", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  skip_if_not_installed("airway")

  se <- airway_se(n_genes = 40)
  design <- stats::model.matrix(
    ~ dex + cell,
    data = as.data.frame(SummarizedExperiment::colData(se))
  )
  shared <- brmDE:::edger_shared_loglik(
    SummarizedExperiment::assay(se, "counts"),
    design
  )
  spline.pts <- log2(shared$phi / 0.1)
  expect_equal(range(spline.pts), c(-20, 10))
  expect_equal(length(shared$phi), 61L)
})

test_that("dispersion_quadratic_log_sd is NA when the curve is not concave", {
  x <- seq(-1, 1, length.out = 11)
  y <- x^2
  expect_true(is.na(
    brmDE:::dispersion_quadratic_log_sd(x, y, x_mode = 0, n_local = 11)
  ))
})

test_that("dispersion_quadratic_log_sd recovers a known Normal SD", {
  sigma <- 0.4
  x_mode <- log(0.1)
  x <- seq(x_mode - 2, x_mode + 2, length.out = 21)
  y <- stats::dnorm(x, mean = x_mode, sd = sigma, log = TRUE)
  expect_equal(
    brmDE:::dispersion_quadratic_log_sd(x, y, x_mode = x_mode, n_local = 5),
    sigma,
    tolerance = 1e-8
  )
})

test_that("dispersion_quadratic_log_sd_over_genes recovers a known SD on a constructed grid", {
  sigma <- 0.35
  n_gene <- 12
  phi <- 0.1 * 2^seq(-10, 10, length.out = 21)
  phi_trend <- withr::with_seed(1, exp(stats::rnorm(n_gene, mean = log(0.1), sd = 0.08)))
  x <- log(phi)
  shared <- t(vapply(seq_len(n_gene), function(g) {
    stats::dnorm(x, mean = log(phi_trend[[g]]), sd = sigma, log = TRUE)
  }, numeric(length(phi))))
  got <- brmDE:::dispersion_quadratic_log_sd_over_genes(
    shared,
    phi,
    prior_n = 1,
    phi_trend = phi_trend,
    n_local = 5
  )
  expect_equal(got, rep(sigma, n_gene), tolerance = 1e-8)
})

# Negative binomial counts with log(phi) ~ N(log(0.1), sigma^2). The Laplace
# SD is the curvature of edgeR's weighted shared likelihood, not sd(log phi)
# itself, but it should track that hyperparameter.
simulate_log_phi_se <- function(sigma, seed, n_gene = 60) {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  skip_if_not_installed("SummarizedExperiment")

  n_rep <- 6L
  group <- factor(rep(c("A", "B"), each = n_rep))
  n_sample <- length(group)
  sim <- withr::with_seed(seed, {
    log_phi <- stats::rnorm(n_gene, mean = log(0.1), sd = sigma)
    mu <- matrix(exp(stats::rnorm(n_gene, 6, 0.3)), n_gene, n_sample)
    mu[, group == "B"] <- mu[, group == "B"] * 1.2
    counts <- matrix(
      stats::rnbinom(n_gene * n_sample, mu = mu, size = 1 / exp(log_phi)),
      n_gene
    )
    list(log_phi = log_phi, counts = counts)
  })
  rownames(sim$counts) <- sprintf("G%03d", seq_len(n_gene))
  colnames(sim$counts) <- sprintf("S%02d", seq_len(n_sample))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = sim$counts),
    colData = data.frame(group = group, row.names = colnames(sim$counts)),
    rowData = data.frame(true_log_phi = sim$log_phi, row.names = rownames(sim$counts))
  )
}

test_that("estimate_dispersion_prior tracks a known log-dispersion SD in simulated counts", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")

  se_tight <- simulate_log_phi_se(sigma = 0.2, seed = 1)
  se_wide <- simulate_log_phi_se(sigma = 0.8, seed = 2)
  se_tight <- estimate_dispersion_prior(se_tight, formula_abundance = ~ group)
  se_wide <- estimate_dispersion_prior(se_wide, formula_abundance = ~ group)

  est_tight <- SummarizedExperiment::rowData(se_tight)$dispersion_prior_log_sd
  est_wide <- SummarizedExperiment::rowData(se_wide)$dispersion_prior_log_sd
  true_tight <- stats::sd(SummarizedExperiment::rowData(se_tight)$true_log_phi)
  true_wide <- stats::sd(SummarizedExperiment::rowData(se_wide)$true_log_phi)
  med_tight <- stats::median(est_tight)
  med_wide <- stats::median(est_wide)

  expect_true(all(is.finite(est_tight)))
  expect_true(all(is.finite(est_wide)))
  expect_lt(abs(med_tight - true_tight) / true_tight, 0.4)
  expect_lt(abs(med_wide - true_wide) / true_wide, 0.4)
  expect_gt(med_wide, med_tight)
})
