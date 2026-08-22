#' Hypothesis tests for one gene fit
#'
#' Thin wrapper around [brms::hypothesis()] that returns a tidy table. Pass
#' arbitrary hypothesis strings for an arbitrary model, or set
#' `hypothesis = "random_vs_rest"` to contrast each random-intercept level
#' against the mean of the others (the tissue test in the immuneBodyMap
#' pipeline).
#'
#' When `group` is set, both population-level (`scope = "standard"`) and
#' group-specific (`scope = "coef"`) tests are returned.
#'
#' @section Point and one-sided hypotheses:
#' The form of the equation decides how `post_prob` and `evid_ratio` are
#' obtained, and the two routes have different requirements.
#'
#' A one-sided hypothesis (`"dextrt > 0"`) is a tally of posterior draws:
#' `post_prob` is the fraction of draws satisfying the inequality, and nothing
#' beyond the posterior is needed.
#'
#' A point hypothesis (`"dextrt = 0"`) is a Bayes factor computed by the
#' Savage-Dickey density ratio, the posterior density at 0 divided by the prior
#' density at 0.
#'
#' That Bayes factor is the Bayesian counterpart of the frequentist
#' likelihood-ratio test, which is what edgeR's `glmLRT()` reports for the same
#' coefficient: both weigh the model with `dextrt` free against the same model
#' with `dextrt` held at 0, so both need the null to be a restriction of the
#' alternative. What differs is the treatment of the parameters. The
#' likelihood-ratio test maximises the likelihood under each hypothesis and
#' refers the ratio to a chi-square distribution to get a p-value; the Bayes
#' factor averages the likelihood over the prior instead, so it needs no null
#' distribution, but it inherits a dependence on that prior. `evid_ratio` is
#' therefore the analogue of the test statistic and `post_prob` its normalised
#' form, on a scale where 1 is no evidence either way. Neither is a p-value:
#' there is no error rate being controlled behind them, so do not threshold
#' them as if there were.
#'
#' Three things follow from how that ratio is computed, and each of them turns
#' `post_prob` and `evid_ratio` into `NA` or into noise when it is not met.
#'
#' * **It needs prior draws.** The fit has to come from [estimate_gene()] with
#'   `sample_prior = "yes"`, which is the default. A fit made with
#'   `sample_prior = "no"` has no prior draws and reports `NA`.
#' * **It needs a proper prior on every parameter tested.** brms cannot draw
#'   from a flat prior, so a parameter left with brms' default flat prior
#'   reports `NA` however many draws the fit has. The default prior set built
#'   by [estimate_gene()] is proper throughout; a `prior` you supply yourself
#'   replaces it entirely, so keep it proper.
#' * **It rests on a density estimate, so it needs draws.** The estimate and
#'   the interval settle long before a density at a single point does. A fit of
#'   a few hundred draws gives a usable `ci_lower`/`ci_upper` next to a
#'   `post_prob` that moves from run to run; around a thousand post-warmup
#'   draws is where it starts to hold still.
#'
#' The ratio also depends on the prior scale rather than on the data alone.
#' Widening `coefficient_prior_scale` in [estimate_gene()] thins prior mass at
#' 0 and so moves the evidence towards the null however clear the data are
#' (Lindley's paradox), which is why that default is one doubling of
#' expression rather than something diffuse.
#'
#' Group-level coefficients are the one case no refit fixes: brms cannot sample
#' their priors, so a point test with `class = "r"` always reports `NA`. Since
#' `hypothesis = "random_vs_rest"` builds `= 0` contrasts of exactly those
#' parameters, read `estimate` and the interval there rather than `post_prob`.
#'
#' @param fit A `brmsfit` from [estimate_gene()].
#' @param hypothesis Character vector of brms hypothesis equations, or the
#'   string `"random_vs_rest"`. An equation with an inequality
#'   (`"dextrt > 0"`) is answered from the posterior alone; one with `=`
#'   (`"dextrt = 0"`) is a point hypothesis and needs prior draws. See
#'   *Point and one-sided hypotheses*.
#' @param group Grouping factor for `scope = "coef"` tests (e.g.
#'   `"tissue_groups"`). Ignored when `NULL`.
#' @param grouping For `hypothesis = "random_vs_rest"`, the random-effect
#'   grouping name (defaults to `group`).
#' @param par Random-effect coefficient to contrast when using
#'   `"random_vs_rest"`. Default `"Intercept"`.
#' @param class Passed to [brms::hypothesis()] (e.g. `"r"` for random
#'   effects).
#' @param robust,alpha Passed to [brms::hypothesis()].
#' @param ... Additional arguments passed to [brms::hypothesis()].
#'
#' @return A tibble with `component`, `group`, `hypothesis`, `estimate`,
#'   interval, and posterior probability columns.
#'
#' @examples
#' \dontrun{
#' data("airway", package = "airway")
#' se <- airway["ENSG00000120129", ]
#' se$offset <- log(colSums(SummarizedExperiment::assay(airway, "counts")))
#' se$dex <- relevel(factor(se$dex), ref = "untrt")
#' fit <- estimate_gene(se, ~ dex + (1 | cell), offset = "offset")
#'
#' # Point hypothesis: needs the prior draws that estimate_gene() takes by
#' # default, and reports a Savage-Dickey Bayes factor in `evid_ratio`.
#' hypothesis_gene(fit, "dextrt = 0")
#'
#' # One-sided: a tally of posterior draws, so no prior draws are involved.
#' hypothesis_gene(fit, "dextrt > 0")
#'
#' # Contrasts of group-level coefficients: `estimate` and the interval are the
#' # answer here, since brms cannot sample priors for these parameters.
#' hypothesis_gene(
#'   fit,
#'   hypothesis = "random_vs_rest",
#'   grouping = "cell",
#'   class = "r"
#' )
#' }
#'
#' @export
hypothesis_gene <- function(fit,
                            hypothesis,
                            group = NULL,
                            grouping = group,
                            par = "Intercept",
                            class = NULL,
                            robust = TRUE,
                            alpha = 0.05,
                            ...) {
  if (missing(hypothesis) || is.null(hypothesis)) {
    stop("`hypothesis` is required.", call. = FALSE)
  }

  vs_rest <- length(hypothesis) == 1L && identical(hypothesis, "random_vs_rest")
  if (vs_rest) {
    if (is.null(grouping)) {
      stop(
        "`grouping` (or `group`) is required when hypothesis = \"random_vs_rest\".",
        call. = FALSE
      )
    }
    hypothesis <- random_intercept_vs_rest_equations(
      fit,
      grouping = grouping,
      par = par
    )
    if (is.null(class)) {
      class <- "r"
    }
    hyp <- brms::hypothesis(
      fit,
      hypothesis,
      class = class,
      robust = robust,
      alpha = alpha,
      ...
    )
    return(tidy_hypothesis(hyp, component = "random", grouping_label = grouping))
  }

  run_one <- function(scope, component, grouping_label, group_arg) {
    args <- list(
      x = fit,
      hypothesis = hypothesis,
      robust = robust,
      alpha = alpha,
      scope = scope
    )
    args$class <- class
    if (!is.null(group_arg)) {
      args$group <- group_arg
    }
    dots <- list(...)
    args[names(dots)] <- dots
    hyp <- do.call(brms::hypothesis, args)
    tidy_hypothesis(hyp, component = component, grouping_label = grouping_label)
  }

  out <- run_one("standard", "fixed", "population", NULL)
  if (!is.null(group) && !vs_rest) {
    out <- dplyr::bind_rows(
      out,
      run_one("coef", "total", group, group)
    )
  }
  out
}
