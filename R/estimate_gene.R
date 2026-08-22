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
#'   `~1` by default. Used only when `shape_prior = "student_t"`, where it
#'   becomes a `shape ~ <terms> + offset(...)` submodel. The offset is
#'   `log(1/<dispersion>)` when `dispersion` is supplied, and `0` when it is
#'   `NULL`, so the terms you give describe departures from that external
#'   estimate or, with a zero offset, the log-shape itself. Do not write the
#'   offset yourself.
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
#'   gene's counts are then free to update. [tidybulk::estimate_dispersion()]
#'   writes this column; this function does not compute it.
#' @param dispersion_degrees_freedom Optional name of the effective degrees
#'   of freedom column that accompanies `dispersion`, as written by
#'   [tidybulk::estimate_dispersion()]. Default `NULL` gives the Student-t
#'   shape intercept a log-scale SD of 1, so the prior scale does not depend
#'   on previous tools. Independent of `dispersion`: you can pass one, both,
#'   or neither.
#' @param shape_prior How the external dispersion is turned into a prior on
#'   the shape parameter. The gene-wise posterior is not constrained to equal
#'   that estimate: both forms only locate and scale the prior, and the
#'   likelihood can move the shape away from it. Both imply the same
#'   log-scale spread, `trigamma(d_eff / 2)`, but they are not
#'   reparameterisations of one another: they differ in tail weight (how
#'   readily a gene may diverge), and they centre different summaries of the
#'   shape on the external estimate. Both put all prior mass on a positive
#'   shape, one by exponentiating an unbounded parameter and the other by
#'   bounding it.
#'
#'   * `"student_t"` (default) adds a
#'     `shape ~ 1 + offset(log(1/dispersion))` submodel (or `offset(0)` when
#'     `dispersion` is `NULL`) and puts a Student-t prior on its intercept,
#'     scaled by `shape_prior_df`. brms gives the submodel a log link, so this
#'     is symmetric in `log(shape)` and centres the *median* shape on
#'     `1/dispersion` (or on `1` when the offset is 0). Heavier-tailed, so
#'     more forgiving when edgeR has over-shrunk a gene, and the submodel
#'     stays available for dispersion covariates.
#'   * `"gamma"` adds no submodel and instead puts
#'     `gamma(d_eff/2, d_eff * dispersion/2)` on `shape` directly. This is the
#'     conjugate form implied by the edgeR hierarchy, in which the dispersion
#'     is scaled inverse chi-square and hence the precision `1/dispersion` is
#'     gamma. It centres the *mean* shape on `1/dispersion`, which sits
#'     `digamma(d_eff/2) - log(d_eff/2)`, roughly `-1/d_eff`, from the
#'     Student-t centre on the log scale. Lighter-tailed. Needs both
#'     `dispersion` and `dispersion_degrees_freedom`, or neither (in which
#'     case brms' vague `gamma(0.01, 0.01)` is used). A scalar `shape`
#'     carries no linear predictor, so a `formula_dispersion` with terms is an
#'     error here.
#' @param shape_prior_df Degrees of freedom \eqn{\nu} of the Student-t prior
#'   on the shape intercept, used only when `shape_prior = "student_t"`. This
#'   single value sets both the `df` argument of `student_t(df, 0, scale)` and
#'   the `scale`, because a Student-t has standard deviation
#'   `scale * sqrt(nu / (nu - 2))`; the scale is therefore
#'   `sd * sqrt((nu - 2) / nu)` for a target `sd`. Must exceed 2, since the
#'   Student-t standard deviation does not exist at or below 2 degrees of
#'   freedom. The target `sd` is `sqrt(trigamma(d_eff / 2))` when
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
#'   This is not only a shrinkage choice. [hypothesis_gene()] evaluates a point
#'   hypothesis such as `"dextrt = 0"` by the Savage-Dickey density ratio, the
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
#'   test a proper prior; brms cannot draw from a flat one, and point
#'   hypotheses on a parameter with a flat prior report NA.
#' @param sample_prior Passed to [brms::brm()]. `"yes"` draws from the prior
#'   alongside the posterior, which is what lets [hypothesis_gene()] report
#'   `pH0` and `evid_ratio` for a point hypothesis such as
#'   `"dextrt = 0"`. The draws cost no sampling work, being taken in generated
#'   quantities, but they do enlarge every stored fit, which at transcriptome
#'   scale is disk you would rather keep. The default `"no"` omits them, since
#'   the default test in [hypothesis_gene()] is directional and needs no prior
#'   draws. Set this to `"yes"` when you intend to ask for a Bayes factor.
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
                          iter = 800,
                          warmup = 300,
                          backend = c("cmdstanr", "rstan"),
                          cores = chains,
                          init = "gene",
                          sanitize_names = FALSE,
                          feature = NULL,
                          stanvars = NULL,
                          ...) {
  backend <- match.arg(backend)
  if (identical(backend, "cmdstanr")) {
    check_and_install_cmdstanr()
  }
  offset <- check_offset_name(offset)
  shape_prior <- match.arg(shape_prior)
  shape_prior_df <- check_student_df(shape_prior_df)
  sample_prior <- match.arg(sample_prior)
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
      dispersion = dispersion,
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
