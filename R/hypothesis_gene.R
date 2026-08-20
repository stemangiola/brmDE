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
