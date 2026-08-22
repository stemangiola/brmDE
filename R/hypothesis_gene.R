#' Hypothesis tests for one gene fit
#'
#' Thin wrapper around [brms::hypothesis()] that returns a tidy table. Pass
#' one or more effects to test; a vector gives one row per entry.
#'
#' When `group` is set, both population-level (`scope = "standard"`) and
#' group-specific (`scope = "coef"`) tests are returned; the second is a
#' [hypothesis_gene_scope()] call stacked under the first.
#'
#' @section One vocabulary for every test:
#' Whichever route a row took, it reports the same two things, so tables from
#' different tests can be read and ranked the same way.
#'
#' `estimate`, the interval, and `log2_fold_change` describe the *contrast*, and
#' nothing else. They do not move when you change the question asked about it,
#' so `"dextrt"` and `"dextrt = 0"` report the same fold change for `dextrt`.
#'
#' `pH0` is the posterior probability that the null being tested is true, and
#' `evid_ratio` is the posterior odds against it, `(1 - pH0) / pH0`. Lower `pH0`
#' and higher `evid_ratio` always mean more evidence for a discovery. What
#' counts as the null differs by route, and `test` records which route a row
#' took, because the routes do not rest on the same assumptions: a `"threshold"`
#' row needs neither a prior nor prior odds, while a `"point"` row needs a
#' proper prior and assumes the null and the alternative were equally likely
#' beforehand. For a point equation, `evid_ratio` is the reciprocal of the Bayes
#' factor [brms::hypothesis()] reports, which measures evidence *for* the null.
#'
#' @section Directional threshold test (the default):
#' Write the effect on its own, `"dextrt"`, and the test asks whether that
#' effect exceeds `test_above_log2FC` in magnitude. Each side is counted
#' separately from the posterior draws, and `pH0` is the probability that the
#' better supported side is wrong:
#'
#' \deqn{p_{H0} = P(\mathrm{ROPE}) + P(\mathrm{losing\ side})
#'   = 1 - \max(p_{+}, p_{-}).}
#'
#' Both ways of being wrong are in there. The effect may be too small to care
#' about, or it may run the other way. Testing `P(|effect| > threshold)`
#' instead would collapse the two and call a posterior sitting on zero with
#' heavy tails a discovery, because such a posterior has mass beyond the
#' threshold on both sides at once while its sign is not even determined. The
#' form above cannot do that: it is at least 1/2 for any posterior symmetric
#' about zero, whatever its tails, and reaches 1 for a posterior lying entirely
#' inside the threshold. At `test_above_log2FC = 0` it reduces to the local
#' false sign rate of Stephens (2017).
#'
#' Which way the effect goes is left to `estimate` and `log2_fold_change`, with
#' `p_positive` and `p_negative` giving the two sides separately. There is no
#' direction label, because a posterior lying entirely inside the threshold
#' leaves both sides at zero and supports no direction at all.
#'
#' Because `pH0` is a probability rather than a test statistic, it turns into
#' a false discovery rate by averaging: see [false_discovery_rate()], which
#' [brmDE()] pipelines apply across genes for you.
#'
#' Nothing here needs prior draws, a proper prior, or a density estimate, so
#' this route has none of the failure modes below and is stable at a few
#' hundred draws.
#'
#' @section Asking for a Bayes factor instead:
#' Write the effect against zero, `"dextrt = 0"`, and `pH0` is filled from a
#' Bayes factor computed by [brms::hypothesis()] rather than counted from the
#' draws. The effect size columns are unaffected, as above.
#'
#' Zero is the only right-hand side either form accepts. Anything else
#' (`"dextrt = 0.7"`, `"dextrt > 0.7"`) is an error rather than a translation,
#' because brms returns draws of `left - right` and a non-zero right-hand side
#' would shift them, reporting a fold change offset by it. Use
#' `test_above_log2FC` to say how large an effect has to be.
#'
#' That Bayes factor is computed by the Savage-Dickey density ratio, the posterior density at 0 divided by the prior
#' density at 0. That Bayes factor is the Bayesian counterpart of the
#' frequentist likelihood-ratio test, which is what edgeR's `glmLRT()` reports
#' for the same coefficient: both weigh the model with `dextrt` free against
#' the same model with `dextrt` held at 0, so both need the null to be a
#' restriction of the alternative. What differs is the treatment of the
#' parameters. The likelihood-ratio test maximises the likelihood under each
#' hypothesis and refers the ratio to a chi-square distribution to get a
#' p-value; the Bayes factor averages the likelihood over the prior instead, so
#' it needs no null distribution, but it inherits a dependence on that prior.
#' `evid_ratio` is the analogue of the test statistic, on a scale where 1 is no
#' evidence either way. It is not a p-value: there is no error rate being
#' controlled behind it, so do not threshold it as if there were.
#'
#' Three things follow from how that ratio is computed, and each of them turns
#' `pH0` and `evid_ratio` into `NA` or into noise when it is not met.
#'
#' * **It needs prior draws.** The fit has to come from [estimate_gene()] with
#'   `sample_prior = "yes"`, which is *not* the default: prior draws enlarge
#'   every stored fit and the directional test above has no use for them, so a
#'   fit made without them reports `NA` here and has to be refitted.
#' * **It needs a proper prior on every parameter tested.** brms cannot draw
#'   from a flat prior, so a parameter left with brms' default flat prior
#'   reports `NA` however many draws the fit has. The default prior set built
#'   by [estimate_gene()] is proper throughout; a `prior` you supply yourself
#'   replaces it entirely, so keep it proper.
#' * **It rests on a density estimate, so it needs draws.** The estimate and
#'   the interval settle long before a density at a single point does. A fit of
#'   a few hundred draws gives a usable `ci_lower`/`ci_upper` next to a `pH0`
#'   that moves from run to run; around a thousand post-warmup draws is where
#'   it starts to hold still.
#'
#' The ratio also depends on the prior scale rather than on the data alone.
#' Widening `coefficient_prior_scale` in [estimate_gene()] thins prior mass at
#' 0 and so moves the evidence towards the null however clear the data are
#' (Lindley's paradox), which is why that default is one doubling of
#' expression rather than something diffuse.
#'
#' Group-level coefficients are the one case no refit fixes: brms cannot sample
#' their priors, so a point equation with `class = "r"` always reports `NA`.
#' The directional test has no such limitation.
#'
#' @references
#' Stephens, M. (2017). False discovery rates: a new deal. *Biostatistics*
#' 18(2), 275-294.
#'
#' @param fit A `brmsfit` from [estimate_gene()].
#' @param hypothesis Character vector of effects to test. An entry written on
#'   its own (`"dextrt"`, or a contrast such as `"dextrt - dexb"`) gets the
#'   directional threshold test. An entry written as `"dextrt = 0"` gets a
#'   Bayes factor instead. The two may be mixed in one call; any other
#'   right-hand side is an error.
#' @param test_above_log2FC Magnitude an effect has to exceed, in log2 fold
#'   change units, so the default `1` asks for a doubling of expression and `2`
#'   for a quadrupling. `0` tests only the sign. Ignored for `"= 0"` entries.
#' @param group Grouping factor for `scope = "coef"` tests (e.g.
#'   `"tissue_groups"`). Ignored when `NULL`.
#' @param class Passed to [brms::hypothesis()] (e.g. `"r"` for random
#'   effects).
#' @param robust,alpha Passed to [brms::hypothesis()].
#' @param ... Additional arguments passed to [brms::hypothesis()].
#'
#' @return A tibble with one row per effect tested and the same columns
#'   whichever route was taken: `component`, `group`, `hypothesis` (the string
#'   you wrote, or its name), `estimate`, `ci_lower`, `ci_upper` on the
#'   natural-log scale of the coefficients, `log2_fold_change` in log2 units,
#'   `test` naming the route (`"threshold"` or `"point"`), `p_positive` and
#'   `p_negative` for the two sides of a threshold test and `NA` otherwise,
#'   `pH0`, and `evid_ratio`. See *One vocabulary for every test*.
#'
#' @examples
#' \dontrun{
#' data("airway", package = "airway")
#' se <- airway["ENSG00000120129", ]
#' se$offset <- log(colSums(SummarizedExperiment::assay(airway, "counts")))
#' se$dex <- relevel(factor(se$dex), ref = "untrt")
#' fit <- estimate_gene(se, ~ dex + (1 | cell), offset = "offset")
#'
#' # Default: is the effect bigger than one log2 fold change, and which way?
#' hypothesis_gene(fit, "dextrt")
#'
#' # Sign only, which is the local false sign rate.
#' hypothesis_gene(fit, "dextrt", test_above_log2FC = 0)
#'
#' # Each cell-line intercept, as four group-level tests.
#' hypothesis_gene(
#'   fit,
#'   c(
#'     "`cell[N052611,Intercept]`",
#'     "`cell[N061011,Intercept]`",
#'     "`cell[N080611,Intercept]`",
#'     "`cell[N61311,Intercept]`"
#'   ),
#'   class = "r"
#' )
#'
#' # Point null: a Savage-Dickey Bayes factor. Needs sample_prior = "yes".
#' hypothesis_gene(fit, "dextrt = 0")
#' }
#'
#' @export
hypothesis_gene <- function(fit,
                            hypothesis,
                            test_above_log2FC = 1,
                            group = NULL,
                            class = NULL,
                            robust = TRUE,
                            alpha = 0.05,
                            ...) {
  out <- hypothesis_gene_scope(
    fit,
    hypothesis,
    test_above_log2FC = test_above_log2FC,
    scope = "standard",
    component = "fixed",
    grouping_label = "population",
    class = class,
    robust = robust,
    alpha = alpha,
    ...
  )
  if (!is.null(group)) {
    out <- dplyr::bind_rows(
      out,
      hypothesis_gene_scope(
        fit,
        hypothesis,
        test_above_log2FC = test_above_log2FC,
        scope = "coef",
        component = "total",
        grouping_label = group,
        group_arg = group,
        class = class,
        robust = robust,
        alpha = alpha,
        ...
      )
    )
  }
  out
}

#' Hypothesis tests at one brms scope
#'
#' Evaluates the same hypotheses either as population-level coefficients
#' (`scope = "standard"`) or as group-specific ones (`scope = "coef"`).
#' [hypothesis_gene()] always runs the population-level scope; when `group` is
#' set it stacks a second call at `scope = "coef"`.
#'
#' @param fit A `brmsfit` from [estimate_gene()].
#' @param hypothesis Effects to test, as in [hypothesis_gene()].
#' @param test_above_log2FC Magnitude an effect has to exceed, as in
#'   [hypothesis_gene()].
#' @param scope Passed to [brms::hypothesis()]: `"standard"` for population-level
#'   coefficients, `"coef"` for group-specific ones.
#' @param component Value written to the returned `component` column.
#' @param grouping_label Value written to the returned `group` column when brms
#'   does not supply one.
#' @param group_arg Grouping factor passed to [brms::hypothesis()] when
#'   `scope = "coef"`. `NULL` otherwise.
#' @param class,robust,alpha,... Passed to [brms::hypothesis()].
#'
#' @return A tibble with the same columns as [hypothesis_gene()].
#'
#' @export
hypothesis_gene_scope <- function(fit,
                                  hypothesis,
                                  test_above_log2FC = 1,
                                  scope = "standard",
                                  component = "fixed",
                                  grouping_label = "population",
                                  group_arg = NULL,
                                  class = NULL,
                                  robust = TRUE,
                                  alpha = 0.05,
                                  ...) {
  point <- hypothesis_is_point(hypothesis)

  # One brms call answers both routes. Whichever of the two supported forms an
  # entry takes, its right-hand side is 0, so the draws brms hands back for
  # `left - right` are the contrast itself: the effect size can be read off the
  # same call that returns the point null's Bayes factor, and a threshold test
  # needs nothing further than those draws.
  hyp <- call_hypothesis(
    fit,
    ifelse(point, hypothesis, paste0("(", hypothesis, ") > 0")),
    robust = robust,
    alpha = alpha,
    class = class,
    scope = scope,
    group_arg = group_arg,
    dots = list(...)
  )
  draws <- hyp$samples

  # scope = "coef" expands one expression into one row per level of the
  # grouping factor, blocked by expression, so anything held per expression
  # repeats in blocks rather than cycling.
  per_expression <- ncol(draws) / length(hypothesis)
  point <- rep(point, each = per_expression)
  labels <- rep(hypothesis_labels(hypothesis), each = per_expression)

  # USE.NAMES would carry brms' internal "H1", "H2" labels into the columns.
  over_draws <- function(f) vapply(draws, f, numeric(1), USE.NAMES = FALSE)
  location <- over_draws(if (isTRUE(robust)) stats::median else mean)
  threshold <- test_above_log2FC * log(2)
  p_positive <- over_draws(function(d) mean(d > threshold))
  p_negative <- over_draws(function(d) mean(d < -threshold))

  # pH0 is the probability of the null on both routes: brms' Post.Prob already
  # is that for a point equation, and the threshold test keeps whichever of the
  # two directional claims the draws support better. evid_ratio is then the
  # posterior odds against the null, so higher always means more evidence for a
  # discovery; for a point equation it is the reciprocal of the Bayes factor
  # brms reports, which measures evidence *for* the null.
  posterior <- if ("Post.Prob" %in% names(hyp$hypothesis)) {
    hyp$hypothesis$Post.Prob
  } else {
    NA_real_
  }
  pH0 <- ifelse(point, posterior, 1 - pmax(p_positive, p_negative))

  tibble::tibble(
    component = component,
    group = if ("Group" %in% names(hyp$hypothesis)) {
      hyp$hypothesis$Group
    } else {
      grouping_label
    },
    hypothesis = labels,
    estimate = location,
    ci_lower = over_draws(function(d) {
      stats::quantile(d, alpha / 2, names = FALSE)
    }),
    ci_upper = over_draws(function(d) {
      stats::quantile(d, 1 - alpha / 2, names = FALSE)
    }),
    log2_fold_change = location / log(2),
    test = ifelse(point, "point", "threshold"),
    p_positive = ifelse(point, NA_real_, p_positive),
    p_negative = ifelse(point, NA_real_, p_negative),
    pH0 = pH0,
    evid_ratio = (1 - pH0) / pH0
  )
}

# Two forms are supported, and both put 0 on the right-hand side of brms'
# `left - right`: a bare contrast, which takes the threshold test, and
# `contrast = 0`, which takes the Bayes factor. That shared zero is what makes
# the returned draws the contrast itself, so the effect size is a property of
# the contrast alone and does not move with the question asked about it. A
# non-zero right-hand side would shift those draws and report a fold change
# offset by it, so it is refused rather than quietly translated.
hypothesis_is_point <- function(hypothesis) {
  at <- regexpr("[=<>]", hypothesis)
  right <- suppressWarnings(as.numeric(trimws(substring(hypothesis, at + 1L))))
  point <- at > 0L & substr(hypothesis, at, at) == "=" &
    !is.na(right) & right == 0
  if (any(at > 0L & !point)) {
    stop(
      "Unsupported hypothesis: ",
      paste0("\"", hypothesis[at > 0L & !point], "\"", collapse = ", "),
      ".\nWrite the contrast on its own (\"dextrt\") to test whether it ",
      "exceeds test_above_log2FC, or \"dextrt = 0\" for a Bayes factor ",
      "against the point null.",
      call. = FALSE
    )
  }
  unname(point)
}

# The label brms would print is built from the parsed equation, which says
# nothing useful once the comparison is ours rather than the user's, so the
# expression the user wrote (or the name they gave it) is used instead.
hypothesis_labels <- function(hypothesis) {
  nms <- names(hypothesis)
  if (is.null(nms)) {
    return(unname(hypothesis))
  }
  ifelse(!is.na(nms) & nzchar(nms), nms, unname(hypothesis))
}

call_hypothesis <- function(fit, hypothesis, robust, alpha, class, scope,
                            group_arg, dots) {
  args <- list(
    x = fit,
    hypothesis = unname(hypothesis),
    robust = robust,
    alpha = alpha,
    scope = scope
  )
  args$class <- class
  if (!is.null(group_arg)) {
    args$group <- group_arg
  }
  args[names(dots)] <- dots
  do.call(brms::hypothesis, args)
}

#' False discovery rate from per-gene null probabilities
#'
#' Turns the `pH0` column of [hypothesis_gene()] into a false discovery rate
#' across genes. Genes are ordered by `pH0` and each gene is given the mean
#' `pH0` of every gene ordered at or before it, so a gene's value is the
#' expected proportion of false discoveries in the set obtained by cutting
#' there. The returned vector is in the order it was given, not in sorted
#' order.
#'
#' This follows from `pH0` being a probability rather than a test statistic.
#' The expected number of false discoveries in a selected set is the sum of the
#' null probabilities of its members, so dividing by the size of the set gives
#' the expected proportion directly, with no null distribution and no p-values
#' in between. The cumulative mean of an ascending sequence is already
#' non-decreasing, so unlike Benjamini-Hochberg q-values no monotonicity step
#' is needed.
#'
#' Two things to keep in mind. The value is a property of the gene set, so
#' filtering genes moves every gene's rate, exactly as it does for
#' Benjamini-Hochberg. And `pH0` is counted from posterior draws, so its
#' resolution is `1 / ndraws`: a `pH0` of 0 means "below what this many draws
#' can resolve", and a rate far below `1 / ndraws` is not supported by the fit
#' it came from.
#'
#' @param pH0 Numeric vector of per-gene null probabilities, as returned in
#'   the `pH0` column of [hypothesis_gene()].
#'
#' @return A numeric vector the same length as `pH0`.
#'
#' @examples
#' pH0 <- c(0.4, 0.001, 0.02, 0.9, 0.005)
#' false_discovery_rate(pH0)
#'
#' @export
false_discovery_rate <- function(pH0) {
  pH0 <- as.numeric(pH0)
  dplyr::cummean(pH0[order(pH0)])[order(order(pH0))]
}

# The gene-wise tables come back one per gene, but a false discovery rate is a
# property of the whole set, so the tables are stacked, the rate computed
# within each hypothesis across genes, and the rows handed back in the shape
# they arrived in.
add_hypothesis_fdr <- function(tables) {
  usable <- length(tables) > 1L &&
    all(vapply(
      tables,
      function(x) is.data.frame(x) && "pH0" %in% names(x),
      logical(1)
    ))
  if (!usable) {
    return(tables)
  }
  sizes <- vapply(tables, nrow, integer(1))
  dplyr::bind_rows(tables) |>
    dplyr::group_by(.data$component, .data$hypothesis) |>
    dplyr::mutate(fdr = false_discovery_rate(.data$pH0)) |>
    dplyr::ungroup() |>
    split_rows(sizes)
}


