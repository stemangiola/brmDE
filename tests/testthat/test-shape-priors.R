default_zinb_priors <- function(data,
                                abundance,
                                offset,
                                formula = NULL,
                                dispersion_degrees_freedom = NULL,
                                shape_prior = "student_t",
                                shape_prior_df = 3) {
  brmDE:::default_gene_priors(
    data = data,
    formula = formula,
    abundance = abundance,
    offset = offset,
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
    formula = f,
    dispersion_degrees_freedom = "dispersion_degrees_freedom",
    shape_prior_df = 10
  )
  shape_int <- priors$prior[priors$class == "Intercept" & priors$dpar == "shape"]
  expect_identical(shape_int, "student_t(10, 0, brmde_shape_scale)")
  expect_equal(
    make_default_zinb_stanvar(
      "brmde_shape_scale",
      dat, "counts", "offset",
      formula = f,
      dispersion_degrees_freedom = "dispersion_degrees_freedom",
      shape_prior_df = 10
    ),
    sqrt(trigamma(9.81176 / 2)) * sqrt(8 / 10),
    tolerance = 1e-6
  )
})

labels_of <- function(formula) attr(stats::terms(formula), "term.labels")

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

test_that("formula_dispersion must be one-sided", {
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
})

test_that("prepare_formula rejects an unknown shape_prior", {
  expect_error(
    brmDE:::prepare_formula(
      ~dex,
      abundance = "counts",
      offset = "offset",
      shape_prior = "lognormal"
    ),
    "should be one of"
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

test_that("shape_gamma_parameters is the conjugate gamma on the deviation", {
  d_eff <- 9.81176
  dat <- tibble::tibble(dispersion_degrees_freedom = d_eff)
  pars <- brmDE:::shape_gamma_parameters(dat, "dispersion_degrees_freedom")
  # The gene enters through the offset, not the prior, so both parameters are
  # d_eff/2 and the deviation exp(Intercept_shape) has prior mean 1.
  expect_equal(pars$shape, d_eff / 2)
  expect_equal(pars$rate, d_eff / 2)
  expect_equal(pars$shape / pars$rate, 1)
  # ... and CV sqrt(2/d), the familiar large-d approximation.
  expect_equal(sqrt(pars$shape) / pars$shape, sqrt(2 / d_eff))
})

test_that("the offset form induces the scalar conjugate gamma on the shape", {
  # shape = exp(Intercept_shape) / phi, and a gamma is closed under scaling,
  # so Gamma(d_eff/2, d_eff/2) on the deviation is Gamma(d_eff/2,
  # d_eff * phi / 2) on the shape: the form the edgeR hierarchy implies.
  for (d_eff in c(4, 9.81176, 30)) {
    for (phi in c(0.01, 0.05, 0.4)) {
      pars <- brmDE:::shape_gamma_parameters(
        tibble::tibble(dispersion_degrees_freedom = d_eff),
        "dispersion_degrees_freedom"
      )
      expect_equal(pars$shape, d_eff / 2)
      expect_equal(pars$rate * phi, d_eff * phi / 2)
      # Mean of the induced prior is edgeR's own point estimate.
      expect_equal(pars$shape / (pars$rate * phi), 1 / phi)
    }
  }
})

test_that("the gamma and Student-t routes imply the same log-scale spread", {
  # chi^2_d is Gamma(d/2, scale = 2), so Var(log X) = trigamma(d/2) either way.
  # Scaling by 1/phi shifts the log scale without changing its spread.
  for (d_eff in c(4, 9.81176, 30)) {
    dat <- tibble::tibble(dispersion_degrees_freedom = d_eff)
    pars <- brmDE:::shape_gamma_parameters(dat, "dispersion_degrees_freedom")
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
      tibble::tibble(dispersion_degrees_freedom = d_eff),
      "dispersion_degrees_freedom"
    )
    gamma_log_centre <- digamma(pars$shape) - log(pars$rate * phi)
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

test_that('shape_prior = "gamma" requires dispersion_degrees_freedom', {
  expect_error(
    brmDE:::shape_gamma_parameters(tibble::tibble(), NULL),
    "requires `dispersion_degrees_freedom`"
  )
  dat <- airway_one_gene_tbl()
  expect_error(
    estimate_gene(dat, ~dex, offset = "offset", shape_prior = "gamma"),
    "requires `dispersion_degrees_freedom`"
  )
})

test_that("shape_gamma_parameters needs a usable degrees of freedom column", {
  expect_error(
    brmDE:::shape_gamma_parameters(
      tibble::tibble(dispersion = 0.05), "dispersion_degrees_freedom"
    ),
    "not found"
  )
  expect_error(
    brmDE:::shape_gamma_parameters(
      tibble::tibble(dispersion_degrees_freedom = Inf),
      "dispersion_degrees_freedom"
    ),
    "finite and positive"
  )
})

test_that("the assembled formula does not depend on shape_prior", {
  args <- list(
    ~dex,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion"
  )
  f_t <- do.call(brmDE:::prepare_formula, c(args, shape_prior = "student_t"))
  f_g <- do.call(brmDE:::prepare_formula, c(args, shape_prior = "gamma"))
  expect_true(brmDE:::has_shape_submodel(f_g))
  expect_equal(
    paste(deparse(f_g$pforms$shape), collapse = " "),
    "shape ~ 1 + offset(log(1/dispersion))"
  )
  expect_equal(deparse(f_t$pforms$shape), deparse(f_g$pforms$shape))
})

test_that('shape_prior = "gamma" puts a flat prior on the shape intercept', {
  d_eff <- 9.81176
  dat <- airway_one_gene_tbl()
  dat$dispersion <- 0.05
  dat$dispersion_degrees_freedom <- d_eff
  f <- brmDE:::prepare_formula(
    ~dex,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion",
    shape_prior = "gamma"
  )
  defaults <- default_zinb_priors(
    dat, "counts", "offset",
    formula = f,
    dispersion_degrees_freedom = "dispersion_degrees_freedom",
    shape_prior = "gamma"
  )
  priors <- defaults$prior
  # Stan has no log-gamma, so the intercept is left flat and the density is
  # added to the target by a stanvar instead. Nothing lands on class "shape",
  # which no longer exists as a scalar parameter.
  shape_int <- priors$prior[priors$class == "Intercept" & priors$dpar == "shape"]
  expect_identical(shape_int, "")
  expect_false(any(priors$class == "shape"))

  expect_equal(stanvar_value(defaults$stanvars, "brmde_shape_gamma_shape"), d_eff / 2)
  expect_equal(stanvar_value(defaults$stanvars, "brmde_shape_gamma_rate"), d_eff / 2)
})

test_that('shape_prior = "gamma" emits the log-gamma target increment', {
  dat <- airway_one_gene_tbl()
  dat$dispersion <- 0.05
  dat$dispersion_degrees_freedom <- 9.81176
  f <- brmDE:::prepare_formula(
    ~dex,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion",
    shape_prior = "gamma"
  )
  scode <- default_zinb_priors(
    dat, "counts", "offset",
    formula = f,
    dispersion_degrees_freedom = "dispersion_degrees_freedom",
    shape_prior = "gamma"
  )$stanvars[["shape_gamma_lprior"]]$scode
  # The gamma is on exp(Intercept_shape); the trailing term is its Jacobian.
  expect_match(scode, "gamma_lpdf\\(exp\\(Intercept_shape\\)", fixed = FALSE)
  expect_match(scode, "\\+ Intercept_shape;")
})

test_that('shape_prior = "gamma" accepts dispersion covariates', {
  dat <- airway_one_gene_tbl()
  dat$dispersion <- 0.05
  dat$dispersion_degrees_freedom <- 9.81176
  f <- brmDE:::prepare_formula(
    ~dex,
    formula_dispersion = ~cell,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion",
    shape_prior = "gamma"
  )
  expect_equal(
    paste(deparse(f$pforms$shape), collapse = " "),
    "shape ~ cell + offset(log(1/dispersion))"
  )
  priors <- make_default_zinb_priors(
    dat, "counts", "offset",
    formula = f,
    dispersion_degrees_freedom = "dispersion_degrees_freedom",
    shape_prior = "gamma"
  )
  # The submodel coefficients need a proper prior of their own; the gamma
  # only speaks to the intercept.
  shape_b <- priors$prior[priors$class == "b" & priors$dpar == "shape"]
  expect_identical(shape_b, "student_t(3, 0, 2)")
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
    formula = f
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
