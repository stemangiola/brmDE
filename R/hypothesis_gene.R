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
#' @param fit A `brmsfit` from [estimate_gene()].
#' @param hypothesis Character vector of brms hypothesis equations, or the
#'   string `"random_vs_rest"`.
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
#'   interval, and posterior probability columns, followed by the convergence
#'   of each contrast.
#'
#' @details
#' Each row carries its own `rhat`, `ess_bulk`, and `mcse`, which
#' [posterior::summarise_draws()] computes from the draws of that contrast
#' rather than from the parameters it was built from ([brms::hypothesis()]
#' returns those draws but no diagnostics for them).
#' The contrast is the quantity being reported, and it commonly mixes better
#' than its parts: a random-effect level that drifts between chains cancels in
#' a difference of levels. `mcse` is the Monte Carlo error of `estimate`, in
#' its units and for the same functional `robust` selects, so it can be read
#' against the width of the interval; `ess_bulk` is the more informative of
#' the two at the few hundred draws a gene-wise fit is usually given.
#'
#' @examples
#' \dontrun{
#' data("airway", package = "airway")
#' se <- airway["ENSG00000120129", ]
#' se$offset <- log(colSums(SummarizedExperiment::assay(airway, "counts")))
#' se$dex <- relevel(factor(se$dex), ref = "untrt")
#' fit <- estimate_gene(se, ~ dex + (1 | cell), offset = "offset")
#'
#' hypothesis_gene(fit, "dextrt = 0")
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

  # Expand the tissue-style shortcut into one equation per random-intercept
  # level, and default to class "r" so brms looks in the random effects.
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
  }

  # Build the brms call by hand so NULL class/group are omitted rather than
  # passed through, and so `...` can still override any named argument.
  run_one <- function(scope, group_arg) {
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
    do.call(brms::hypothesis, args)
  }

  # Population-level tests always run. With `group`, also test the same
  # equations at the group-specific (fixed + random) coefficients.
  hyps <- list(run_one("standard", NULL))
  components <- if (vs_rest) "random" else "fixed"
  grouping_labels <- if (vs_rest) grouping else "population"
  if (!is.null(group) && !vs_rest) {
    hyps <- c(hyps, list(run_one("coef", group)))
    components <- c(components, "total")
    grouping_labels <- c(grouping_labels, group)
  }

  # One tidy pass: rename brms columns, then attach contrast-level rhat,
  # ess_bulk and mcse. The chain count is all that flattening of hyp$samples
  # loses, and the fit is still here.
  n_chains <- brms::nchains(fit)
  dplyr::bind_rows(Map(function(hyp, component, grouping_label) {
    tbl <- tibble::as_tibble(hyp$hypothesis)
    if (!nrow(tbl)) {
      return(tbl)
    }
    out <- dplyr::transmute(
      tbl,
      component = component,
      group = if ("Group" %in% names(tbl)) .data$Group else grouping_label,
      hypothesis = .data$Hypothesis,
      estimate = .data$Estimate,
      ci_lower = .data$CI.Lower,
      ci_upper = .data$CI.Upper,
      post_prob = if ("Post.Prob" %in% names(tbl)) .data$Post.Prob else NA_real_,
      evid_ratio = if ("Evid.Ratio" %in% names(tbl)) .data$Evid.Ratio else NA_real_
    )
    dplyr::bind_cols(out, contrast_diagnostics(hyp, n_chains, robust = robust))
  }, hyps, components, grouping_labels))
}
