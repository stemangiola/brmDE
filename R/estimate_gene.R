#' Fit a brms count model for one gene
#'
#' Prepares one gene's samples (data frame or one-row
#' `SummarizedExperiment`), injects a precomputed offset into the formula if
#' needed, and fits a brms model. Priors, inits, and MCMC defaults follow the
#' immuneBodyMap ZINB specification and can all be overridden.
#'
#' @param data A data frame with one row per sample, or a
#'   `SummarizedExperiment` with a single feature.
#' @param formula A `formula` or `brmsformula`. The response may be omitted
#'   (`~ dex + (1 | donor)`); it is then set to `abundance`. An
#'   `offset(<offset>)` term is added if missing.
#' @param offset Required name of the precomputed offset column in `data`
#'   (or in `colData` if `data` is a `SummarizedExperiment`). The offset is
#'   never calculated inside this function.
#' @param dispersion Optional name of a precomputed edgeR dispersion column
#'   in `rowData` (or in `data`), used to build the prior on the negative
#'   binomial shape as described under `shape_prior`. [estimate_dispersion()]
#'   writes this column; this function does not compute it.
#' @param dispersion_degrees_freedom Name of the effective degrees of freedom
#'   column that accompanies `dispersion`, as written by
#'   [estimate_dispersion()]. Only consulted when `dispersion` is supplied;
#'   when the column is absent the shape prior falls back to a weakly
#'   informative default.
#' @param shape_prior How edgeR's dispersion is turned into a prior on the
#'   shape parameter. Both forms imply the same log-scale spread,
#'   `trigamma(d_eff / 2)`, but they are not reparameterisations of one
#'   another: they differ in tail weight, and they centre different summaries
#'   of the shape on edgeR's estimate. Both put all prior mass on a positive
#'   shape, one by exponentiating an unbounded parameter and the other by
#'   bounding it. See [estimate_dispersion()] for the derivation.
#'
#'   * `"student_t"` (default) adds a
#'     `shape ~ 1 + offset(log(1/dispersion))` submodel and puts a Student-t
#'     prior on its intercept, scaled by `shape_prior_df`. brms gives the
#'     submodel a log link, so this is symmetric in `log(shape)` and centres
#'     the *median* shape on `1/dispersion`. Heavier-tailed, so more forgiving
#'     when edgeR has over-shrunk a gene, and the submodel stays available for
#'     dispersion covariates.
#'   * `"gamma"` adds no submodel and instead puts
#'     `gamma(d_eff/2, d_eff * dispersion/2)` on `shape` directly. This is the
#'     conjugate form implied by the edgeR hierarchy, in which the dispersion
#'     is scaled inverse chi-square and hence the precision `1/dispersion` is
#'     gamma. It centres the *mean* shape on `1/dispersion`, which sits
#'     `digamma(d_eff/2) - log(d_eff/2)`, roughly `-1/d_eff`, from the
#'     Student-t centre on the log scale. Lighter-tailed. Ignored with a
#'     warning if `formula` already carries a shape submodel.
#' @param shape_prior_df Degrees of freedom \eqn{\nu} of the Student-t prior
#'   on the shape intercept, used only when `shape_prior = "student_t"`. This
#'   single value sets both the `df` argument of `student_t(df, 0, scale)` and
#'   the `scale`, because a Student-t has standard deviation
#'   `scale * sqrt(nu / (nu - 2))`; the scale is therefore
#'   `sd * sqrt((nu - 2) / nu)` for a target `sd`. Must exceed 2, since the
#'   Student-t standard deviation does not exist at or below 2 degrees of
#'   freedom. The target `sd` comes from the effective degrees of freedom
#'   recorded by [estimate_dispersion()]; see its Details for the derivation.
#' @param family A brms family. Default is
#'   [brms::zero_inflated_negbinomial()].
#' @param abundance Name of the count column, or assay name if `data` is a
#'   `SummarizedExperiment`.
#' @param prior Optional brms prior. If `NULL` and `family` is ZINB, student-t
#'   priors are used with the intercept centred at
#'   `mean(log1p(counts / exp(offset)))`.
#' @param chains,iter,warmup MCMC settings. Defaults match the manuscript
#'   pipeline (2 chains, 600 iterations, 400 warmup).
#' @param backend Passed to [brms::brm()]. Default `"cmdstanr"`.
#' @param cores Number of CPU cores. Defaults to `min(available cores, chains)`.
#' @param init Inits for [brms::brm()]. `"gene"` (default) builds intercept
#'   and coefficient inits from `make_standata()`. Use `"random"` or a list
#'   for brms defaults / custom inits.
#' @param sanitize_names If `TRUE`, collapse repeated underscores in column
#'   names (`___` → `_`) to work around a brms parsing issue.
#' @param feature Optional gene id used only in a warning string (handy for
#'   targets logs).
#' @param ... Additional arguments passed to [brms::brm()].
#'
#' @return A `brmsfit` object.
#'
#' @examples
#' print("cmdstanr is needed to run this example.")
#' # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
#'
#' if (instantiate::stan_cmdstan_exists()) {
#'   data("airway", package = "airway")
#'   se <- airway["ENSG00000120129", ]
#'   se$offset <- log(colSums(SummarizedExperiment::assay(airway, "counts")))
#'   se$dex <- relevel(factor(se$dex), ref = "untrt")
#'
#'   fit <- estimate_gene(
#'     se,
#'     formula = ~ dex + (1 | cell),
#'     offset = "offset",
#'     abundance = "counts"
#'   )
#' }
#'
#' @export
estimate_gene <- function(data,
                          formula,
                          offset,
                          dispersion = NULL,
                          dispersion_degrees_freedom = "dispersion_degrees_freedom",
                          shape_prior = c("student_t", "gamma"),
                          shape_prior_df = 3,
                          family = NULL,
                          abundance = "counts",
                          prior = NULL,
                          chains = 2,
                          iter = 600,
                          warmup = 400,
                          backend = "cmdstanr",
                          cores = NULL,
                          init = "gene",
                          sanitize_names = FALSE,
                          feature = NULL,
                          ...) {
  if (identical(backend, "cmdstanr")) {
    check_and_install_cmdstanr()
  }
  offset <- check_offset_name(offset)
  shape_prior <- match.arg(shape_prior)
  shape_prior_df <- check_student_df(shape_prior_df)
  if (!is.null(dispersion)) {
    dispersion <- check_dispersion_name(dispersion)
    dispersion_degrees_freedom <-
      check_degrees_freedom_name(dispersion_degrees_freedom)
  } else {
    dispersion_degrees_freedom <- NULL
  }
  prepared <- prepare_gene_data(
    data,
    abundance = abundance,
    offset = offset,
    dispersion = dispersion,
    dispersion_degrees_freedom = dispersion_degrees_freedom,
    sanitize_names = sanitize_names
  )
  data <- prepared$data
  abundance <- prepared$abundance
  offset <- prepared$offset
  dispersion <- prepared$dispersion
  dispersion_degrees_freedom <- prepared$dispersion_degrees_freedom

  if (is.null(family)) {
    family <- brms::zero_inflated_negbinomial()
  }

  if (!is.null(feature)) {
    message(glue::glue("***** Gene:___{feature}___*****"))
  } else if (".feature" %in% names(data)) {
    feature <- unique(data[[".feature"]])
    if (length(feature) == 1L && !is.na(feature)) {
      message(glue::glue("***** Gene:___{feature}___*****"))
    }
  }

  formula <- prepare_formula(
    formula,
    abundance = abundance,
    offset = offset,
    dispersion = dispersion,
    shape_prior = shape_prior
  )

  if (is.null(prior) && is_zinb_family(family)) {

    # Prior for the shape parameter
    if (has_shape_submodel(formula)) {
      prior <- shape_student_t_prior(
        data,
        formula,
        dispersion,
        dispersion_degrees_freedom,
        shape_prior_df,
        shape_prior
      )
    } else if (identical(shape_prior, "gamma")) {
      prior <- shape_gamma_prior(data, dispersion, dispersion_degrees_freedom)
    } else {
      prior <- brms::prior(student_t(3, 0, 2), class = shape)
    }

    # Prior for the location parameters
    prior <- c(prior, zinb_location_priors(data, abundance, offset))
  }

  if (identical(init, "gene")) {
    init <- make_gene_inits(
      formula = formula,
      data = data,
      family = family,
      prior = prior,
      chains = chains,
      abundance = abundance,
      offset = offset
    )
  }

  n_cores <- if (is.null(cores)) min(detect_cores(), chains) else cores
  tpc <- floor(n_cores / chains)
  threads <- if (tpc <= 1L) NULL else brms::threading(tpc)

  res <- brms::brm(
    formula = formula,
    data = data,
    family = family,
    prior = prior,
    chains = chains,
    cores = n_cores,
    threads = threads,
    warmup = warmup,
    iter = iter,
    backend = backend,
    init = init,
    ...
  )

  if (nrow(res$data) != nrow(data)) {
    stop(
      "brms fitted ", nrow(res$data), " of the ", nrow(data),
      " rows supplied; the rest were dropped, most likely by missing values ",
      "in a model term. The fit would not describe the data you passed.",
      call. = FALSE
    )
  }

  res
}
