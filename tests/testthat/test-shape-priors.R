default_zinb_priors <- function(data,
                                abundance,
                                offset,
                                formula = NULL,
                                dispersion_prior_log_sd = NULL,
                                shape_prior_df = 3) {
  brmDE:::default_gene_priors(
    data = data,
    formula = formula,
    abundance = abundance,
    offset = offset,
    dispersion_prior_log_sd = dispersion_prior_log_sd,
    shape_prior_df = shape_prior_df
  )
}

make_default_zinb_priors <- function(...) default_zinb_priors(...)$prior

# Prior constants that depend on the gene reach Stan as data, so their values
# are asserted on the stanvars rather than on the prior string.
stanvar_value <- function(stanvars, name) stanvars[[name]]$sdata

make_default_zinb_stanvar <- function(name, ...) {
  stanvar_value(default_zinb_priors(...)$stanvars, name)
}

test_that("dispersion_log_sd_from_degrees_freedom is the trigamma SD, not the 2/d approximation", {
  # Var(log s^2) = trigamma(d/2) for s^2 ~ sigma^2 chi^2_d / d.
  for (d in c(4, 6, 10, 30)) {
    expect_equal(brmDE:::dispersion_log_sd_from_degrees_freedom(d), sqrt(trigamma(d / 2)))
  }
  # sqrt(2/d) understates it, most severely at low df.
  expect_gt(brmDE:::dispersion_log_sd_from_degrees_freedom(4), sqrt(2 / 4))
  expect_lt(brmDE:::dispersion_log_sd_from_degrees_freedom(150) - sqrt(2 / 150), 0.01)
  # robust = TRUE can shrink completely (prior.df = Inf); a failed design is NA.
  expect_true(is.na(brmDE:::dispersion_log_sd_from_degrees_freedom(Inf)))
  expect_true(is.na(brmDE:::dispersion_log_sd_from_degrees_freedom(NA_real_)))
  expect_true(is.na(brmDE:::dispersion_log_sd_from_degrees_freedom(0)))
  expect_true(is.na(brmDE:::dispersion_log_sd_from_degrees_freedom(-1)))
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

test_that("check_log_sd_value reads a precomputed SD, defaulting only when unasked", {
  d_eff <- 9.81176
  dat <- tibble::tibble(
    dispersion = 0.05,
    dispersion_df_prior_sd = brmDE:::dispersion_log_sd_from_degrees_freedom(d_eff)
  )
  expect_equal(
    brmDE:::check_log_sd_value(dat, "dispersion_df_prior_sd"),
    sqrt(trigamma(d_eff / 2))
  )
  expect_error(
    brmDE:::check_log_sd_value(
      tibble::tibble(dispersion = 0.05),
      "dispersion_df_prior_sd"
    ),
    "not found"
  )
})

test_that("check_log_sd_value rejects unusable log-SDs", {
  for (bad in c(Inf, NA_real_, 0, -1)) {
    expect_error(
      brmDE:::check_log_sd_value(
        tibble::tibble(dispersion = 0.05, dispersion_df_prior_sd = bad),
        "dispersion_df_prior_sd"
      ),
      "finite and positive"
    )
  }
})

test_that("default_zinb_priors sets the shape intercept from d_eff", {
  dat <- airway_one_gene_tbl()
  dat$dispersion <- 0.05
  dat$dispersion_df_prior_sd <- brmDE:::dispersion_log_sd_from_degrees_freedom(9.81176)
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
    dispersion_prior_log_sd = "dispersion_df_prior_sd",
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
    dispersion_prior_log_sd = "dispersion_df_prior_sd",
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
  dat$dispersion_df_prior_sd <- brmDE:::dispersion_log_sd_from_degrees_freedom(9.81176)
  f <- brmDE:::prepare_formula(
    ~ dex,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion"
  )
  priors <- make_default_zinb_priors(
    dat, "counts", "offset",
    formula = f,
    dispersion_prior_log_sd = "dispersion_df_prior_sd",
    shape_prior_df = 10
  )
  shape_int <- priors$prior[priors$class == "Intercept" & priors$dpar == "shape"]
  expect_identical(shape_int, "student_t(10, 0, brmde_shape_scale)")
  expect_equal(
    make_default_zinb_stanvar(
      "brmde_shape_scale",
      dat, "counts", "offset",
      formula = f,
      dispersion_prior_log_sd = "dispersion_df_prior_sd",
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

test_that("a dispersion column without a width column uses the default scale", {
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

test_that("dispersion_quadratic_log_sd recovers a Normal SD from curvature", {
  sigma <- 0.4
  x <- seq(-1, 1, length.out = 11)
  y <- -(x^2) / (2 * sigma^2)
  expect_equal(
    brmDE:::dispersion_quadratic_log_sd(x, y, x_mode = 0, n_local = 11),
    sigma,
    tolerance = 1e-8
  )
})

test_that("curvature uses dispersion_log_sd as the Student-t target SD", {
  dat <- airway_one_gene_tbl()
  dat$dispersion <- 0.05
  dat$dispersion_prior_sd <- 0.35
  f <- brmDE:::prepare_formula(
    ~dex,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion"
  )
  priors <- make_default_zinb_priors(
    dat, "counts", "offset",
    formula = f,
    dispersion_prior_log_sd = "dispersion_prior_sd"
  )
  shape_int <- priors$prior[priors$class == "Intercept" & priors$dpar == "shape"]
  expect_identical(shape_int, "student_t(3, 0, brmde_shape_scale)")
  expect_equal(
    stanvar_value(default_zinb_priors(
      dat, "counts", "offset",
      formula = f,
      dispersion_prior_log_sd = "dispersion_prior_sd"
    )$stanvars, "brmde_shape_scale"),
    brmDE:::student_t_scale_for_sd(0.35, 3)
  )
})

test_that("curvature defaults to log-scale SD 1", {
  dat <- airway_one_gene_tbl()
  f <- suppressMessages(
    brmDE:::prepare_formula(
      ~dex,
      abundance = "counts",
      offset = "offset"
    )
  )
  expect_equal(
    stanvar_value(default_zinb_priors(
      dat, "counts", "offset",
      formula = f
    )$stanvars, "brmde_shape_scale"),
    brmDE:::student_t_scale_for_sd(1, 3)
  )
})

test_that("shape covariates still get a Student-t prior on class b", {
  dat <- airway_one_gene_tbl()
  dat$dispersion <- 0.05
  f <- brmDE:::prepare_formula(
    ~dex,
    formula_dispersion = ~cell,
    abundance = "counts",
    offset = "offset",
    dispersion = "dispersion"
  )
  priors <- make_default_zinb_priors(dat, "counts", "offset", formula = f)
  shape_b <- priors$prior[priors$class == "b" & priors$dpar == "shape"]
  expect_identical(shape_b, "student_t(3, 0, 2)")
})

test_that("as_gene_tibble copies a rowData dispersion column", {
  se <- airway_one_gene()
  SummarizedExperiment::rowData(se)$dispersion <- 0.4
  tbl <- brmDE:::as_gene_tibble(se, abundance = "counts", rowdata_cols = "dispersion")
  expect_equal(unique(tbl$dispersion), 0.4)
  expect_equal(nrow(tbl), ncol(se))
})
