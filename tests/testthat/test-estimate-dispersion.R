test_that("estimate_dispersion writes a positive tagwise column on airway", {
  skip_if_not_installed("edgeR")

  se <- airway_se(n_genes = 150)
  out <- estimate_dispersion(se, ~ dex + (1 | cell), abundance = "counts")
  disp <- SummarizedExperiment::rowData(out)$dispersion
  expect_equal(length(disp), nrow(se))
  expect_true(all(is.finite(disp)))
  expect_true(all(disp > 0))
})

default_zinb_priors <- function(data,
                                abundance,
                                offset,
                                formula = NULL,
                                dispersion = NULL,
                                dispersion_degrees_freedom = NULL,
                                shape_prior = "student_t",
                                shape_prior_df = 3) {
  brmDE:::default_gene_priors(
    data = data,
    formula = formula,
    abundance = abundance,
    offset = offset,
    dispersion = dispersion,
    dispersion_degrees_freedom = dispersion_degrees_freedom,
    shape_prior_df = shape_prior_df,
    shape_prior = shape_prior
  )
}

make_default_zinb_priors <- function(...) default_zinb_priors(...)$prior

# Prior constants that depend on the gene reach Stan as data, so their values
# are asserted on the stanvars rather than on the prior string.
stanvar_value <- function(stanvars, name) stanvars[[name]]$sdata

make_default_zinb_stanvar <- function(name, ...) {
  stanvar_value(default_zinb_priors(...)$stanvars, name)
}

test_that("estimate_dispersion records d_eff = df.residual + prior.df", {
  skip_if_not_installed("edgeR")

  se <- airway_se(n_genes = 150)
  out <- estimate_dispersion(se, ~ dex + (1 | cell), abundance = "counts")
  d_eff <- SummarizedExperiment::rowData(out)$dispersion_degrees_freedom
  expect_equal(length(d_eff), nrow(se))

  # (1 | cell) enters the edgeR design as a fixed `cell`.
  design <- stats::model.matrix(
    ~ dex + cell,
    data = droplevels(as.data.frame(SummarizedExperiment::colData(se)))
  )
  fit <- edgeR::estimateDisp(
    SummarizedExperiment::assay(se, "counts"),
    design = design
  )
  expect_equal(
    unique(d_eff),
    (ncol(se) - ncol(design)) + fit$prior.df
  )
})

test_that("dispersion_log_sd is the trigamma SD, not the 2/d approximation", {
  # Var(log s^2) = trigamma(d/2) for s^2 ~ sigma^2 chi^2_d / d.
  for (d in c(4, 6, 10, 30)) {
    expect_equal(brmDE:::dispersion_log_sd(d), sqrt(trigamma(d / 2)))
  }
  # sqrt(2/d) understates it, most severely at low df.
  expect_gt(brmDE:::dispersion_log_sd(4), sqrt(2 / 4))
  expect_lt(brmDE:::dispersion_log_sd(150) - sqrt(2 / 150), 0.01)
})

test_that("student_t_scale_for_sd inverts the Student-t SD formula", {
  for (nu in c(3, 5, 10)) {
    sd_target <- 0.475
    scale <- brmDE:::student_t_scale_for_sd(sd_target, nu)
    # Student-t SD is scale * sqrt(nu / (nu - 2)).
    expect_equal(scale * sqrt(nu / (nu - 2)), sd_target)
  }
  expect_equal(brmDE:::student_t_scale_for_sd(1, 3), 1 / sqrt(3))
})

test_that("check_student_df requires df above 2", {
  expect_equal(brmDE:::check_student_df(3), 3)
  expect_error(brmDE:::check_student_df(2), "greater than 2")
  expect_error(brmDE:::check_student_df(1), "greater than 2")
  expect_error(brmDE:::check_student_df(Inf), "greater than 2")
  expect_error(brmDE:::check_student_df(c(3, 4)), "greater than 2")
})

test_that("shape_intercept_scale uses d_eff, defaulting only when unasked", {
  d_eff <- 9.81176
  dat <- tibble::tibble(dispersion = 0.05, dispersion_degrees_freedom = d_eff)
  expect_equal(
    brmDE:::shape_intercept_scale(dat, "dispersion_degrees_freedom", nu = 3),
    sqrt(trigamma(d_eff / 2)) / sqrt(3)
  )

  # Degrees of freedom omitted: documented default SD, converted to scale.
  expect_equal(
    brmDE:::shape_intercept_scale(dat, NULL, nu = 3),
    brmDE:::student_t_scale_for_sd(brmDE:::shape_prior_sd_default, 3)
  )
  # A df name that is not a column is an error, not a silent default.
  expect_error(
    brmDE:::shape_intercept_scale(
      tibble::tibble(dispersion = 0.05), "dispersion_degrees_freedom", nu = 3
    ),
    "not found"
  )
})

test_that("shape_intercept_scale rejects unusable degrees of freedom", {
  # robust = TRUE can shrink completely, giving prior.df = Inf, and a design
  # edgeR could not fit gives NA. Neither is quietly replaced by a default.
  for (bad in c(Inf, NA_real_, 0, -1)) {
    expect_error(
      brmDE:::shape_intercept_scale(
        tibble::tibble(dispersion = 0.05, dispersion_degrees_freedom = bad),
        "dispersion_degrees_freedom",
        nu = 3
      ),
      "finite and positive"
    )
  }
})

test_that("estimate_dispersion output column names are arguments", {
  skip_if_not_installed("edgeR")

  se <- airway_se(n_genes = 150)
  out <- estimate_dispersion(
    se,
    ~ dex + (1 | cell),
    abundance = "counts",
    dispersion = "phi",
    dispersion_degrees_freedom = "phi_deff"
  )
  rd <- SummarizedExperiment::rowData(out)
  expect_true(all(c("phi", "phi_deff") %in% names(rd)))
  expect_false(any(c("dispersion", "dispersion_degrees_freedom") %in% names(rd)))
  expect_true(all(rd$phi > 0))
  expect_length(unique(rd$phi_deff), 1L)

  expect_error(
    estimate_dispersion(se, ~dex, dispersion = "x", dispersion_degrees_freedom = "x"),
    "different columns"
  )
})

test_that("default_zinb_priors sets the shape intercept from d_eff", {
  dat <- airway_one_gene_tbl()
  dat$dispersion <- 0.05
  dat$dispersion_degrees_freedom <- 9.81176
  f <- brmDE:::prepare_formula(
    ~ dex,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion"
  )
  priors <- make_default_zinb_priors(
    dat,
    "counts",
    "offset",
    formula = f,
    dispersion = "dispersion",
    dispersion_degrees_freedom = "dispersion_degrees_freedom",
    shape_prior_df = 3
  )
  shape_int <- priors$prior[priors$class == "Intercept" & priors$dpar == "shape"]
  expect_length(shape_int, 1L)
  expect_identical(shape_int, "student_t(3, 0, brmde_shape_scale)")

  # The scale itself travels as data, so that the Stan code is the same for
  # every gene and one compiled model can serve all of them.
  scale <- make_default_zinb_stanvar(
    "brmde_shape_scale",
    dat,
    "counts",
    "offset",
    formula = f,
    dispersion = "dispersion",
    dispersion_degrees_freedom = "dispersion_degrees_freedom",
    shape_prior_df = 3
  )
  expect_equal(scale, sqrt(trigamma(9.81176 / 2)) * sqrt(1 / 3), tolerance = 1e-6)

  # shape ~ 1 + offset(...) has no population-level coefficients, so a
  # `b` prior on shape would not correspond to any model parameter.
  expect_false(any(priors$class == "b" & priors$dpar == "shape"))
})

test_that("shape_prior_df changes both the df and the scale", {
  dat <- airway_one_gene_tbl()
  dat$dispersion <- 0.05
  dat$dispersion_degrees_freedom <- 9.81176
  f <- brmDE:::prepare_formula(
    ~ dex,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion"
  )
  priors <- make_default_zinb_priors(
    dat, "counts", "offset",
    formula = f, dispersion = "dispersion",
    dispersion_degrees_freedom = "dispersion_degrees_freedom",
    shape_prior_df = 10
  )
  shape_int <- priors$prior[priors$class == "Intercept" & priors$dpar == "shape"]
  expect_identical(shape_int, "student_t(10, 0, brmde_shape_scale)")
  expect_equal(
    make_default_zinb_stanvar(
      "brmde_shape_scale",
      dat, "counts", "offset",
      formula = f, dispersion = "dispersion",
      dispersion_degrees_freedom = "dispersion_degrees_freedom",
      shape_prior_df = 10
    ),
    sqrt(trigamma(9.81176 / 2)) * sqrt(8 / 10),
    tolerance = 1e-6
  )
})

labels_of <- function(formula) attr(stats::terms(formula), "term.labels")

test_that("fixed_effects_formula turns random intercepts into fixed factors", {
  # 1:group simplifies to group.
  expect_equal(
    labels_of(brmDE:::fixed_effects_formula(~ dex + (1 | cell))),
    c("dex", "cell")
  )
  expect_equal(
    labels_of(brmDE:::fixed_effects_formula(~ (1 | cell))),
    "cell"
  )
  # An intercept-only formula still has nothing to put in the design.
  expect_equal(labels_of(brmDE:::fixed_effects_formula(~1)), character(0))
})

test_that("fixed_effects_formula expands random slopes against the group", {
  # (1 + f1:f2 | cell) is 1:cell + f1:f2:cell, i.e. cell + f1:f2:cell.
  # terms() then orders each interaction by where its variables first appear,
  # so the label comes back as cell:f1:f2.
  expect_equal(
    labels_of(brmDE:::fixed_effects_formula(~ dex + (1 + f1:f2 | cell))),
    c("dex", "cell", "cell:f1:f2")
  )
  expect_equal(
    labels_of(brmDE:::fixed_effects_formula(~ (dex | cell))),
    c("cell", "cell:dex")
  )
  # No random intercept, so no bare grouping factor.
  expect_equal(
    labels_of(brmDE:::fixed_effects_formula(~ (0 + dex | cell))),
    "dex:cell"
  )
  # `||` is the uncorrelated form of the same fixed structure.
  expect_equal(
    labels_of(brmDE:::fixed_effects_formula(~ (dex || cell))),
    c("cell", "cell:dex")
  )
  # brms' (terms | ID | group) only ties posteriors together across
  # distributional parameters; the fixed structure is that of (terms | group).
  expect_equal(
    labels_of(brmDE:::fixed_effects_formula(~ dex + (1 | p | cell))),
    c("dex", "cell")
  )
  expect_equal(
    labels_of(brmDE:::fixed_effects_formula(~ (dex | p | cell))),
    c("cell", "cell:dex")
  )
  # Nested grouping expands the way a formula would.
  expect_equal(
    labels_of(brmDE:::fixed_effects_formula(~ (1 | donor / cell))),
    c("donor", "donor:cell")
  )
})

test_that("estimate_dispersion fails when the expanded design is saturated", {
  skip_if_not_installed("edgeR")

  se <- simulated_se()
  # One level per sample, so `(1 | sid)` as a fixed effect uses up the design.
  se$sid <- factor(colnames(se))
  expect_error(
    estimate_dispersion(se, ~ treatment + (1 | sid)),
    "no residual degrees of freedom"
  )
  expect_error(
    estimate_dispersion(se[, 0], ~treatment),
    "no samples"
  )
})

test_that("estimate_dispersion expands a random slope on a larger dataset", {
  skip_if_not_installed("edgeR")

  se <- simulated_se()
  out <- estimate_dispersion(se, ~ treatment + (treatment | donor))
  disp <- SummarizedExperiment::rowData(out)$dispersion
  d_eff <- SummarizedExperiment::rowData(out)$dispersion_degrees_freedom
  expect_true(all(is.finite(disp) & disp > 0))

  # (treatment | donor) becomes donor + treatment:donor, so the design costs
  # 12 of the 24 samples rather than the 2 a fixed-effects-only formula would.
  design <- stats::model.matrix(
    ~ treatment + donor + treatment:donor,
    data = as.data.frame(SummarizedExperiment::colData(se))
  )
  expect_equal(ncol(design), 12L)
  fit <- edgeR::estimateDisp(
    SummarizedExperiment::assay(se, "counts"),
    design = design
  )
  expect_equal(unique(d_eff), (ncol(se) - ncol(design)) + fit$prior.df)

  # A random intercept alone spends fewer degrees of freedom than a slope.
  d_eff_intercept <- SummarizedExperiment::rowData(
    estimate_dispersion(se, ~ treatment + (1 | donor))
  )$dispersion_degrees_freedom
  expect_gt(unique(d_eff_intercept), unique(d_eff))
})

test_that("estimate_dispersion recovers the dispersion it was simulated with", {
  skip_if_not_installed("edgeR")

  se <- simulated_se()
  out <- estimate_dispersion(se, ~ treatment + (1 | donor))
  expect_gt(
    stats::cor(
      SummarizedExperiment::rowData(out)$dispersion,
      SummarizedExperiment::rowData(out)$true_dispersion,
      method = "spearman"
    ),
    0.6
  )
})

test_that("fixed_effects_formula does not duplicate a term already fixed", {
  expect_equal(
    labels_of(brmDE:::fixed_effects_formula(~ cell + (1 | cell))),
    "cell"
  )
})

test_that("estimate_dispersion rejects non-SE input", {
  expect_error(estimate_dispersion(mtcars, ~ dex), "SummarizedExperiment")
})

test_that("prepare_formula adds a shape offset from dispersion", {
  f <- brmDE:::prepare_formula(
    ~ dex,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion"
  )
  txt <- paste(deparse(f$formula), collapse = " ")
  shape_txt <- paste(deparse(f$pforms$shape), collapse = " ")
  expect_match(txt, "offset\\s*\\(\\s*offset\\s*\\)")
  expect_match(shape_txt, "log\\s*\\(\\s*1/dispersion")
})

test_that("prepare_formula rejects a shape submodel inside formula_abundance", {
  # Silently keeping it would drop the edgeR dispersion offset.
  f0 <- brms::bf(counts ~ dex + offset(offset), shape ~ 1)
  expect_error(
    brmDE:::prepare_formula(
      f0,
      abundance = "counts",
      offset = "offset",
      dispersion = "dispersion"
    ),
    "Model the dispersion through `formula_dispersion`"
  )
})

test_that("formula_dispersion becomes the shape submodel, offset and all", {
  f <- brmDE:::prepare_formula(
    ~dex,
    formula_dispersion = ~cell,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion"
  )
  shape_txt <- paste(deparse(f$pforms$shape), collapse = " ")
  expect_match(shape_txt, "^shape ~ cell \\+ offset\\(log\\(1/dispersion\\)\\)$")

  # Dispersion covariates are fitted even without an edgeR estimate; the
  # offset is then 0 rather than log(1/dispersion).
  f_no_disp <- brmDE:::prepare_formula(
    ~dex,
    formula_dispersion = ~cell,
    abundance = "counts",
    offset = "offset"
  )
  expect_equal(
    paste(deparse(f_no_disp$pforms$shape), collapse = " "),
    "shape ~ cell + offset(0)"
  )
})

test_that("formula_dispersion must be one-sided and gamma-compatible", {
  expect_error(
    brmDE:::prepare_formula(
      ~dex,
      formula_dispersion = shape ~ cell,
      abundance = "counts",
      offset = "offset",
      dispersion = "dispersion"
    ),
    "must be one-sided"
  )
  expect_error(
    brmDE:::prepare_formula(
      ~dex,
      formula_dispersion = ~cell,
      abundance = "counts",
      offset = "offset",
      dispersion = "dispersion",
      shape_prior = "gamma"
    ),
    "scalar `shape` with no linear predictor"
  )
})

test_that("prepare_formula reports the formulas it assembled", {
  expect_message(
    brmDE:::prepare_formula(
      ~dex,
      abundance = "counts",
      offset = "offset",
      dispersion = "dispersion"
    ),
    "Abundance model \\(offset added by brmDE\\): counts ~ dex \\+ offset\\(offset\\)"
  )
  expect_message(
    brmDE:::prepare_formula(
      ~dex,
      abundance = "counts",
      offset = "offset",
      dispersion = "dispersion"
    ),
    "Dispersion model \\(offset added by brmDE\\): shape ~ 1 \\+ offset\\(log\\(1/dispersion\\)\\)"
  )
  expect_message(
    brmDE:::prepare_formula(
      ~dex,
      abundance = "counts",
      offset = "offset"
    ),
    "Dispersion model \\(offset added by brmDE\\): shape ~ 1 \\+ offset\\(0\\)"
  )
  # Nothing is announced as added when the user wrote the offset themselves.
  expect_message(
    brmDE:::prepare_formula(
      counts ~ dex + offset(offset),
      abundance = "counts",
      offset = "offset"
    ),
    "Abundance model: counts ~ dex \\+ offset\\(offset\\)"
  )
})

test_that("shape_gamma_parameters is the conjugate gamma with mean 1/phi", {
  d_eff <- 9.81176
  phi <- 0.05
  dat <- tibble::tibble(dispersion = phi, dispersion_degrees_freedom = d_eff)
  pars <- brmDE:::shape_gamma_parameters(dat, "dispersion", "dispersion_degrees_freedom")
  expect_equal(pars$shape, d_eff / 2)
  expect_equal(pars$rate, d_eff * phi / 2)
  # Gamma(a, rate = b) has mean a/b, here edgeR's own point estimate.
  expect_equal(pars$shape / pars$rate, 1 / phi)
  # ... and CV sqrt(2/d), the familiar large-d approximation.
  expect_equal(sqrt(pars$shape) / pars$shape, sqrt(2 / d_eff))
})

test_that("the gamma and Student-t routes imply the same log-scale spread", {
  # chi^2_d is Gamma(d/2, scale = 2), so Var(log X) = trigamma(d/2) either way.
  for (d_eff in c(4, 9.81176, 30)) {
    dat <- tibble::tibble(dispersion = 0.05, dispersion_degrees_freedom = d_eff)
    pars <- brmDE:::shape_gamma_parameters(dat, "dispersion", "dispersion_degrees_freedom")
    expect_equal(trigamma(pars$shape), brmDE:::dispersion_log_sd(d_eff)^2)
  }
})

test_that("the two routes centre different summaries of the shape", {
  # The Student-t is symmetric in log(shape) under brms' log link, so it
  # centres the median on 1/phi; the gamma centres the mean. The log-scale
  # centres therefore sit digamma(d/2) - log(d/2) apart, roughly -1/d.
  phi <- 0.05
  for (d_eff in c(4, 9.81176, 30)) {
    pars <- brmDE:::shape_gamma_parameters(
      tibble::tibble(dispersion = phi, dispersion_degrees_freedom = d_eff),
      "dispersion",
      "dispersion_degrees_freedom"
    )
    gamma_log_centre <- digamma(pars$shape) - log(pars$rate)
    student_log_centre <- log(1 / phi)
    expect_equal(
      gamma_log_centre - student_log_centre,
      digamma(d_eff / 2) - log(d_eff / 2)
    )
    # The gamma centre always sits below the Student-t one.
    expect_lt(gamma_log_centre, student_log_centre)
  }

  # digamma(d/2) - log(d/2) expands as -1/d - 1/(3 d^2). The leading -1/d
  # alone understates the gap by 8% at d = 4, in the same way sqrt(2/d)
  # understates the spread; two terms are within 0.2% throughout.
  gap <- function(d) digamma(d / 2) - log(d / 2)
  for (d_eff in c(4, 9.81176, 30, 100)) {
    expect_equal(gap(d_eff), -1 / d_eff - 1 / (3 * d_eff^2), tolerance = 2e-3)
  }
  expect_gt(abs(gap(4) - (-1 / 4)) / abs(gap(4)), 0.07)
})

test_that("prepare_formula uses a zero shape offset when dispersion is omitted", {
  f <- suppressMessages(
    brmDE:::prepare_formula(
      ~dex,
      abundance = "counts",
      offset = "offset"
    )
  )
  expect_true(brmDE:::has_shape_submodel(f))
  expect_equal(
    paste(deparse(f$pforms$shape), collapse = " "),
    "shape ~ 1 + offset(0)"
  )
})

test_that("shape_gamma_parameters uses brms' vague gamma when unasked", {
  vague <- list(shape = 0.01, rate = 0.01)
  dat <- tibble::tibble(dispersion = 0.05, dispersion_degrees_freedom = 9.81176)
  expect_equal(brmDE:::shape_gamma_parameters(dat, NULL, NULL), vague)
})

test_that("shape_gamma_parameters needs both columns or neither", {
  dat <- tibble::tibble(dispersion = 0.05, dispersion_degrees_freedom = 9.8)
  expect_error(
    brmDE:::shape_gamma_parameters(dat, "dispersion", NULL),
    "both `dispersion` and `dispersion_degrees_freedom`"
  )
  expect_error(
    brmDE:::shape_gamma_parameters(dat, NULL, "dispersion_degrees_freedom"),
    "both `dispersion` and `dispersion_degrees_freedom`"
  )
  expect_error(
    brmDE:::shape_gamma_parameters(
      tibble::tibble(dispersion = 0.05), "dispersion", "dispersion_degrees_freedom"
    ),
    "not found"
  )
  expect_error(
    brmDE:::shape_gamma_parameters(
      tibble::tibble(dispersion_degrees_freedom = 9.8), "dispersion",
      "dispersion_degrees_freedom"
    ),
    "was not found"
  )
})

test_that("shape_gamma_parameters rejects unusable dispersion or df", {
  expect_error(
    brmDE:::shape_gamma_parameters(
      tibble::tibble(dispersion = 0.05, dispersion_degrees_freedom = Inf),
      "dispersion",
      "dispersion_degrees_freedom"
    ),
    "finite and positive"
  )
  expect_error(
    brmDE:::shape_gamma_parameters(
      tibble::tibble(dispersion = NA_real_, dispersion_degrees_freedom = 9.8),
      "dispersion",
      "dispersion_degrees_freedom"
    ),
    "non-finite or non-positive"
  )
})

test_that('shape_prior = "gamma" builds no shape submodel', {
  f <- brmDE:::prepare_formula(
    ~ dex,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion",
    shape_prior = "gamma"
  )
  expect_false(brmDE:::has_shape_submodel(f))
  expect_match(
    paste(deparse(if (inherits(f, "brmsformula")) f$formula else f), collapse = " "),
    "offset\\s*\\(\\s*offset\\s*\\)"
  )
})

test_that('shape_prior = "gamma" puts a gamma on shape directly', {
  d_eff <- 9.81176
  phi <- 0.05
  dat <- airway_one_gene_tbl()
  dat$dispersion <- phi
  dat$dispersion_degrees_freedom <- d_eff
  f <- brmDE:::prepare_formula(
    ~ dex,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion",
    shape_prior = "gamma"
  )
  priors <- make_default_zinb_priors(
    dat, "counts", "offset",
    formula = f, dispersion = "dispersion",
    dispersion_degrees_freedom = "dispersion_degrees_freedom",
    shape_prior = "gamma"
  )
  shape_prior <- priors$prior[priors$class == "shape"]
  expect_length(shape_prior, 1L)
  expect_identical(
    shape_prior,
    "gamma(brmde_shape_gamma_shape, brmde_shape_gamma_rate)"
  )
  stanvars <- default_zinb_priors(
    dat, "counts", "offset",
    formula = f, dispersion = "dispersion",
    dispersion_degrees_freedom = "dispersion_degrees_freedom",
    shape_prior = "gamma"
  )$stanvars
  expect_equal(stanvar_value(stanvars, "brmde_shape_gamma_shape"), d_eff / 2)
  expect_equal(stanvar_value(stanvars, "brmde_shape_gamma_rate"), d_eff * phi / 2)
  expect_false(any(priors$dpar == "shape"))
})

test_that('shape_prior = "gamma" rejects a user shape submodel', {
  dat <- airway_one_gene_tbl()
  dat$dispersion <- 0.05
  dat$dispersion_degrees_freedom <- 9.81176
  f <- brms::bf(counts ~ dex + offset(offset), shape ~ 1)
  expect_error(
    make_default_zinb_priors(
      dat, "counts", "offset",
      formula = f, dispersion = "dispersion",
      dispersion_degrees_freedom = "dispersion_degrees_freedom",
      shape_prior = "gamma"
    ),
    "has a shape submodel"
  )
})

test_that("a dispersion column without degrees of freedom uses the default scale", {
  dat <- airway_one_gene_tbl()
  dat$dispersion <- 0.05
  f <- brmDE:::prepare_formula(
    ~dex,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion"
  )
  scale <- make_default_zinb_stanvar(
    "brmde_shape_scale",
    dat, "counts", "offset",
    formula = f, dispersion = "dispersion"
  )
  expect_equal(
    scale,
    brmDE:::student_t_scale_for_sd(brmDE:::shape_prior_sd_default, 3)
  )
})

test_that("omitting dispersion still builds a shape intercept prior", {
  dat <- airway_one_gene_tbl()
  f <- suppressMessages(
    brmDE:::prepare_formula(~dex, abundance = "counts", offset = "offset")
  )
  priors <- make_default_zinb_priors(dat, "counts", "offset", formula = f)
  shape_int <- priors$prior[priors$class == "Intercept" & priors$dpar == "shape"]
  expect_identical(shape_int, "student_t(3, 0, brmde_shape_scale)")
  expect_equal(
    make_default_zinb_stanvar(
      "brmde_shape_scale",
      dat, "counts", "offset",
      formula = f
    ),
    brmDE:::student_t_scale_for_sd(brmDE:::shape_prior_sd_default, 3)
  )
})

test_that("check_shape_prior accepts only the two parameterisations", {
  expect_equal(brmDE:::check_shape_prior("student_t"), "student_t")
  expect_equal(brmDE:::check_shape_prior("gamma"), "gamma")
  expect_equal(brmDE:::check_shape_prior(c("student_t", "gamma")), "student_t")
  expect_error(brmDE:::check_shape_prior("normal"))
  expect_error(brmDE:::check_shape_prior(3), "student_t")
})

test_that("as_gene_tibble copies a rowData dispersion column", {
  se <- airway_one_gene()
  SummarizedExperiment::rowData(se)$dispersion <- 0.4
  tbl <- brmDE:::as_gene_tibble(se, abundance = "counts", rowdata_cols = "dispersion")
  expect_equal(unique(tbl$dispersion), 0.4)
  expect_equal(nrow(tbl), ncol(se))
})
