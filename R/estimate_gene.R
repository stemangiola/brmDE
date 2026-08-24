#' Fit a brms count model for one gene
#'
#' Prepares one gene's samples (data frame or one-row
#' `SummarizedExperiment`) and fits a brms model. You give the mean and
#' dispersion models as plain one-sided formulas; the library size offset and the externally estimated dispersion offset are
#' appended here, and the assembled formulas are reported with a message. A
#' precomputed dispersion (typically from [tidybulk::estimate_dispersion()])
#' is used as a **prior** on the negative
#' binomial shape, not as a fixed plug-in: the gene-wise likelihood can still
#' pull the posterior away from that estimate. Priors, inits, and MCMC
#' defaults follow the immuneBodyMap ZINB specification and can all be
#' overridden.
#'
#' Formula terms are resolved against `data` first and the global environment
#' after. The environment a formula was written in is deliberately dropped:
#' R attaches the defining frame to every formula, so a fit produced inside a
#' function would otherwise serialise that whole frame - the
#' `SummarizedExperiment` included - into each stored `brmsfit`, once per gene.
#' A formula referring to a local variable (`~ poly(dose, k)` with a local `k`)
#' will therefore not resolve it; promote such values to arguments or to the
#' global environment. The targets pipeline has always behaved this way, because
#' formulas make a round trip through disk as text.
#'
#' @param data A data frame with one row per sample, or a
#'   `SummarizedExperiment` with a single feature.
#' @param formula_abundance Model for the mean. A `formula` or `brmsformula`
#'   whose response may be omitted (`~ dex + (1 | donor)`); it is then set to
#'   `abundance`. The library size `offset(<offset>)` term is added here rather
#'   than by you, and the assembled formula is reported with a message.
#' @param formula_dispersion One-sided model for the negative binomial shape,
#'   `~1` by default. It becomes a `shape ~ <terms> + offset(...)` submodel
#'   under either `shape_prior`. The offset is `log(1/<dispersion>)` when
#'   `dispersion` is supplied, and `0` when it is `NULL`, so the terms you
#'   give describe departures from that external estimate or, with a zero
#'   offset, the log-shape itself. Do not write the offset yourself.
#' @param offset Required name of the precomputed offset column in `data`
#'   (or in `colData` if `data` is a `SummarizedExperiment`). The offset is
#'   never calculated inside this function.
#' @param dispersion Optional name of a precomputed dispersion column
#'   in `rowData` (or in `data`). Default `NULL` puts a zero offset on the
#'   shape submodel: the intercept is then `log(shape)`, with prior median
#'   shape `1`. Use this when the sample size is large enough that an
#'   informative location from previous tools is unnecessary. When supplied,
#'   the column becomes `offset(log(1/<dispersion>))` as described under
#'   `shape_prior`: that is the prior *location* of the shape, which the
#'   gene's counts are then free to update. Pass
#'   [tidybulk::estimate_dispersion()]'s `dispersion_trended` column: that is
#'   the across-gene trend \eqn{s_0^2}, information this gene has not yet
#'   contributed. Do not pass `dispersion_shrinked`, the tagwise posterior
#'   \eqn{q_g^{post}}, which already includes this gene's counts. See
#'   `vignette("dispersion-priors")`.
#' @param dispersion_degrees_freedom Optional name of the degrees of freedom
#'   column that accompanies `dispersion`, as written by
#'   [tidybulk::estimate_dispersion()] (`dispersion_degrees_freedom`).
#'   Default `NULL` targets a log-scale SD of 1 under `"student_t"`.
#'   `"gamma"` requires this column: it reads \(d_0/2\) as its parameters and
#'   has no fallback. Independent of `dispersion`: under `"student_t"` you
#'   can pass one, both, or neither.
#' @param shape_prior How the external dispersion is turned into a prior on
#'   the shape parameter. The gene-wise posterior is not constrained to equal
#'   that estimate: both forms only locate and scale the prior, and the
#'   likelihood can move the shape away from it.
#'
#'   Both use the same `shape ~ <formula_dispersion> + offset(...)` submodel,
#'   so `dispersion` enters as an offset on brms' log link either way and the
#'   assembled formula does not depend on this argument. Both therefore put
#'   all prior mass on a positive shape by exponentiating, and both accept
#'   dispersion covariates. Both also imply the same log-scale spread,
#'   `trigamma(d_0 / 2)`. What differs is the prior on the submodel
#'   intercept, and hence tail weight and which summary of the shape is
#'   centred on the external estimate. They are not reparameterisations of
#'   one another.
#'
#'   * `"student_t"` (default) puts a Student-t prior on the intercept,
#'     scaled by `shape_prior_df`. That is symmetric in `log(shape)`, so it
#'     centres the *median* shape on `1/dispersion` (or on `1` when the
#'     offset is 0). Heavier-tailed, so more forgiving when edgeR has
#'     over-shrunk a gene.
#'   * `"gamma"` requires `dispersion_degrees_freedom`. It puts
#'     `gamma(d_0/2, d_0/2)` on `exp(intercept)`, the
#'     shape's multiplicative deviation from the offset. Because a gamma is
#'     closed under scaling, the induced prior on the shape itself is exactly
#'     `gamma(d_0/2, d_0 * dispersion/2)`: the conjugate form implied by
#'     the edgeR hierarchy, in which the dispersion is scaled inverse
#'     chi-square and hence the precision `1/dispersion` is gamma. It centres
#'     the *mean* shape on `1/dispersion`, which sits
#'     `digamma(d_0/2) - log(d_0/2)`, roughly `-1/d_0`, from the
#'     Student-t centre on the log scale. Lighter-tailed.
#'
#'   Stan has no log-gamma density, so the gamma form cannot be written as a
#'   `brms` prior on the intercept. It is added to the target directly through
#'   a `stanvar`, which means [brms::prior_summary()] reports the intercept as
#'   flat; the gamma is nonetheless the only prior acting on it, and its two
#'   parameters are recoverable from `fit$stanvars`.
#'
#'   See `vignette("dispersion-priors")` for the scaled-chi-square
#'   hierarchy, why the prior location is the trend \eqn{s_0^2}
#'   (`dispersion_trended`) rather than the tagwise posterior, the trigamma
#'   identity that sets the shared log-scale spread, and why the two forms
#'   are not reparameterisations of each other.
#' @param shape_prior_df Degrees of freedom \eqn{\nu} of the Student-t prior
#'   on the shape intercept, used only when `shape_prior = "student_t"`. This
#'   single value sets both the `df` argument of `student_t(df, 0, scale)` and
#'   the `scale`, because a Student-t has standard deviation
#'   `scale * sqrt(nu / (nu - 2))`; the scale is therefore
#'   `sd * sqrt((nu - 2) / nu)` for a target `sd`. Must exceed 2, since the
#'   Student-t standard deviation does not exist at or below 2 degrees of
#'   freedom. The target `sd` is `sqrt(trigamma(d_0 / 2))` when
#'   `dispersion_degrees_freedom` is supplied, and 1 when it is omitted.
#' @param coefficient_prior_scale Scale of the Student-t prior on the
#'   population-level coefficients of the abundance model (`class = "b"`).
#'   Those coefficients sit behind a log link, so each one is a natural-log
#'   fold change. The default `0.7` is one log2 fold change, a doubling of
#'   expression, which is the order of the effects a differential expression
#'   analysis looks for; `1.4` is two log2 fold changes, and so on. Unlike
#'   `shape_prior_df`, this is the Stan scale itself rather than a target
#'   standard deviation, so it reaches the prior untouched.
#'
#'   This is not only a shrinkage choice. [hypothesis_gene()] evaluates an
#'   entry such as `"dextrt = 0"` by the Savage-Dickey density ratio, the
#'   posterior density at 0 over the prior density at 0, so a wider prior thins
#'   prior mass at 0 and moves the evidence towards the null however clear the
#'   data are. Keep this on the order of the effects you are looking for.
#' @param coefficient_prior_df Degrees of freedom of that Student-t prior. The
#'   default 3 keeps the tails heavy enough that a large fold change is still
#'   cheap.
#' @param family A brms family. Default is
#'   [brms::zero_inflated_negbinomial()].
#' @param abundance Name of the count column, or assay name if `data` is a
#'   `SummarizedExperiment`.
#' @param prior Optional brms prior. If `NULL` and `family` is a negative
#'   binomial (with or without zero inflation), student-t priors are used with
#'   the intercept centred at `mean(log1p(counts / exp(offset)))` and the
#'   coefficients scaled by `coefficient_prior_scale`. Supplying your own
#'   `prior` replaces that set entirely, so give every parameter you intend to
#'   test a proper prior; brms cannot draw from a flat one, and an `"= 0"`
#'   entry on a parameter with a flat prior reports NA.
#' @param sample_prior Passed to [brms::brm()]. `"yes"` draws from the prior
#'   alongside the posterior, which is what lets [hypothesis_gene()] report
#'   `pH0` and `evid_ratio` for an entry such as
#'   `"dextrt = 0"`. The draws cost no sampling work, being taken in generated
#'   quantities, but they do enlarge every stored fit, which at transcriptome
#'   scale is disk you would rather keep. The default `"no"` omits them, since
#'   the default test in [hypothesis_gene()] is directional and needs no prior
#'   draws. Set this to `"yes"` when you intend to ask for a Bayes factor.
#' @param chains,draws_warmup,draws_sampling MCMC settings, counted in draws
#'   per chain rather than in brms' `iter` (which counts warmup and sampling
#'   together): brms is given `warmup = draws_warmup` and
#'   `iter = draws_warmup + draws_sampling`. The fit therefore keeps
#'   `chains * draws_sampling` draws, and that product is what bounds the
#'   resolution of a `pH0` counted from them, `1 / ndraws`. Defaults are the
#'   manuscript pipeline's: 2 chains of 300 warmup and 500 sampling draws, so
#'   1000 draws kept.
#' @param backend Passed to [brms::brm()]. Default `"cmdstanr"`.
#' @param cores Number of CPU cores. Defaults to `min(available cores, chains)`.
#' @param init Inits for [brms::brm()]. `"gene"` (default) builds intercept
#'   and coefficient inits from `make_standata()`. Use `"random"` or a list
#'   for brms defaults / custom inits.
#' @param sanitize_names If `TRUE`, collapse repeated underscores in column
#'   names (`___` → `_`) to work around a brms parsing issue.
#' @param feature Optional gene id used only in a warning string (handy for
#'   targets logs).
#' @param stanvars Optional [brms::stanvar()] object appended to the ones
#'   brmDE builds for its own priors.
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
#'     formula_abundance = ~ dex + (1 | cell),
#'     offset = "offset",
#'     abundance = "counts"
#'   )
#' }
#'
#' @export
estimate_gene <- function(data,
                          formula_abundance,
                          formula_dispersion = ~1,
                          offset,
                          dispersion = NULL,
                          dispersion_degrees_freedom = NULL,
                          shape_prior = c("student_t", "gamma"),
                          shape_prior_df = 3,
                          coefficient_prior_scale = 0.7, # ~1 log2 fold change
                          coefficient_prior_df = 3,
                          family = NULL,
                          abundance = "counts",
                          prior = NULL,
                          sample_prior = c("no", "yes", "only"),
                          chains = 2,
                          draws_warmup = 300,
                          draws_sampling = 500,
                          backend = c("cmdstanr", "rstan"),
                          cores = chains,
                          init = "gene",
                          sanitize_names = FALSE,
                          feature = NULL,
                          stanvars = NULL,
                          ...) {
  backend <- match.arg(backend)
  offset <- check_offset_name(offset)
  shape_prior <- match.arg(shape_prior)
  shape_prior_df <- check_student_df(shape_prior_df)
  sample_prior <- match.arg(sample_prior)
  
  if (identical(shape_prior, "gamma") && is.null(dispersion_degrees_freedom)) {
    stop(
      'brmDE says: `shape_prior = "gamma"` requires `dispersion_degrees_freedom`. ',
      'Use the default `shape_prior = "student_t"` if you do not have it.',
      call. = FALSE
    )
  }
  if (identical(backend, "cmdstanr")) {
    check_and_install_cmdstanr()
  }
  if (!is.null(dispersion)) {
    dispersion <- check_dispersion_name(dispersion)
  }
  if (!is.null(dispersion_degrees_freedom)) {
    dispersion_degrees_freedom <-
      check_degrees_freedom_name(dispersion_degrees_freedom)
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
    formula_abundance,
    formula_dispersion = formula_dispersion,
    abundance = abundance,
    offset = offset,
    dispersion = dispersion,
    shape_prior = shape_prior
  )

  if (is.null(prior) && is_negbinomial_family(family)) {
    defaults <- default_gene_priors(
      data = data,
      formula = formula,
      abundance = abundance,
      offset = offset,
      dispersion_degrees_freedom = dispersion_degrees_freedom,
      shape_prior_df = shape_prior_df,
      shape_prior = shape_prior,
      coefficient_prior_scale = coefficient_prior_scale,
      coefficient_prior_df = coefficient_prior_df
    )
    prior <- defaults$prior
    stanvars <- combine_stanvars(stanvars, defaults$stanvars)
  }

  if (identical(init, "gene")) {
    init <- make_gene_inits(
      formula = formula,
      data = data,
      family = family,
      prior = prior,
      chains = chains,
      abundance = abundance,
      offset = offset,
      stanvars = stanvars
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
    stanvars = stanvars,
    sample_prior = sample_prior,
    chains = chains,
    cores = n_cores,
    threads = threads,
    warmup = draws_warmup,
    iter = draws_warmup + draws_sampling,
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
