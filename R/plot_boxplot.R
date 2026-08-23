#' Plot one gene across a factor, with a posterior predictive check
#'
#' Boxplot of library-size-normalised counts for a single [estimate_gene()]
#' fit, stratified by a factor in the abundance model. Blue boxplots are
#' posterior predictive draws: if the model is descriptively adequate they
#' should roughly overlay the black boxplots of the observed (or adjusted)
#' counts. The construction follows `sccomp_boxplot()`.
#'
#' The library size always comes out of the plotted counts, through the model
#' rather than by dividing: [adjust_gene()] predicts with the offset zeroed,
#' which for a TMM offset `log(1 / multiplier)` leaves the TMM-normalised
#' abundance. Setting `remove_unwanted_effects = TRUE` drops the random
#' effects and any other covariate as well, keeping `factor` alone, so a dex
#' effect is not mixed with cell-line intercepts.
#'
#' @param fit A `brmsfit` from [estimate_gene()].
#' @param factor Character string naming the factor of interest in the
#'   abundance model, used to stratify the boxplot.
#' @param remove_unwanted_effects If `TRUE`, keep `factor` alone: the random
#'   effects and every other covariate come out of both the counts and the
#'   posterior predictive. Defaults to `FALSE`, where only the library size
#'   is removed.
#' @param number_of_draws Number of posterior predictive draws to overlay
#'   (blue boxplots). Defaults to `100`.
#' @param offset Name of the offset column in `fit$data`.
#' @param feature Optional gene id used as the plot subtitle. The fit does
#'   not carry this: pass it from the pipeline result (`.feature`) or from
#'   the gene you filtered on.
#' @param seed Optional seed passed to [brms::posterior_predict()].
#'
#' @return A `ggplot` object.
#'
#' @seealso [estimate_gene()], [adjust_gene()], [hypothesis_gene()]
#'
#' @examples
#' print("cmdstanr is needed to run this example.")
#' # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
#'
#' \donttest{
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
#'
#'   plot_boxplot(fit, factor = "dex", feature = "ENSG00000120129")
#'   plot_boxplot(fit, factor = "dex", feature = "ENSG00000120129",
#'             remove_unwanted_effects = TRUE)
#' }
#' }
#'
#' @export
plot_boxplot <- function(fit,
                      factor,
                      remove_unwanted_effects = FALSE,
                      number_of_draws = 100,
                      offset = "offset",
                      feature = NULL,
                      seed = NULL) {
  abundance <- gene_response_name(fit)
  data_observed <- gene_plot_observed(
    fit,
    factor = factor,
    abundance = abundance,
    offset = offset,
    remove_unwanted_effects = remove_unwanted_effects
  )

  my_boxplot <- ggplot2::ggplot(data_observed)
  posterior_predictive_check <- gene_plot_posterior_predictive_check(
    fit,
    factor = if (isTRUE(remove_unwanted_effects)) factor else NULL,
    offset = offset,
    number_of_draws = number_of_draws,
    seed = seed
  )
  posterior_predictive_check[[factor]] <-
    as.character(posterior_predictive_check[[factor]])

  my_boxplot <- my_boxplot +
    ggplot2::stat_summary(
      ggplot2::aes(
        x = .data[[factor]],
        y = .data$generated,
        group = .data[[factor]]
      ),
      fun.data = calc_boxplot_stat,
      geom = "boxplot",
      linewidth = 0.2,
      colour = "blue",
      data = posterior_predictive_check
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(
        x = .data[[factor]],
        y = .data$normalised,
        group = .data[[factor]]
      ),
      outlier.shape = NA,
      outlier.colour = NA,
      outlier.size = 0,
      median.linewidth = 0.5,
      linewidth = 0.5,
      fill = NA
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(
        x = .data[[factor]],
        y = .data$normalised,
        shape = .data$is_zero
      ),
      position = ggplot2::position_jitter(height = 0, width = 0.2),
      size = 1.5
    ) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21)) +
    ggplot2::scale_y_continuous(
      transform = S_sqrt_trans(),
      labels = drop_leading_zero
    ) +
    ggplot2::xlab(factor) +
    ggplot2::ylab(
      if (isTRUE(remove_unwanted_effects)) {
        "Adjusted counts"
      } else {
        "Normalised counts"
      }
    ) +
    ggplot2::labs(shape = "Zero count") +
    ggplot2::ggtitle(
      sprintf(
        "Grouped by %s (for multi-factor models, associations could be hardly observable with unidimensional data stratification)",
        factor
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

  if (!is.null(feature) && length(feature) == 1L && !is.na(feature) &&
    nzchar(as.character(feature))) {
    my_boxplot <- my_boxplot + ggplot2::labs(subtitle = as.character(feature))
  }
  my_boxplot
}

#' Response column of a gene fit
#'
#' Reads `fit$formula$resp`, falling back to the left-hand side of the
#' abundance formula when that is empty.
#'
#' @param fit A `brmsfit` from [estimate_gene()].
#'
#' @return Character string naming the count column.
#'
#' @noRd
gene_response_name <- function(fit) {
  resp <- fit$formula$resp
  if (is.null(resp) || !nzchar(resp)) {
    resp <- all.vars(fit$formula$formula)[[1]]
  }
  resp
}

#' Fixed-effect covariates other than the plotted factor
#'
#' Names of abundance-model fixed effects that `adjust_gene()` and the
#' posterior predictive should drop, so the plotted effect is not mixed
#' with the rest.
#'
#' @param fit A `brmsfit` from [estimate_gene()].
#' @param factor Character string naming the factor of interest.
#' @param abundance Name of the count column.
#' @param offset Name of the offset column.
#'
#' @return Character vector of covariate names, possibly empty.
#'
#' @noRd
gene_nuisance_fixed_effects <- function(fit, factor, abundance, offset) {
  bt <- brms::brmsterms(fit$formula)
  fe <- bt$dpars$mu$fe$formula
  if (is.null(fe)) {
    return(character())
  }
  setdiff(all.vars(fe), c(abundance, offset, factor))
}

#' Counts for [plot_boxplot()], with the library size removed
#'
#' [adjust_gene()] values, always predicted with the offset zeroed, so the
#' library size comes out of the counts through the model rather than by
#' dividing them. `remove_unwanted_effects` decides what else goes with it:
#' `FALSE` keeps every other effect, `TRUE` keeps `factor` alone. Adds
#' `normalised` and `is_zero`.
#'
#' @param fit A `brmsfit` from [estimate_gene()].
#' @param factor Character string naming the factor of interest.
#' @param abundance Name of the count column.
#' @param offset Name of the offset column.
#' @param remove_unwanted_effects If `TRUE`, also drop the group-level
#'   effects and every covariate other than `factor`.
#'
#' @return The model data with `normalised` and `is_zero` columns.
#'
#' @noRd
gene_plot_observed <- function(fit,
                               factor,
                               abundance,
                               offset,
                               remove_unwanted_effects) {
  keep_factor_only <- isTRUE(remove_unwanted_effects)
  adjusted <- adjust_gene(
    fit,
    nullify = if (keep_factor_only) {
      gene_nuisance_fixed_effects(fit, factor, abundance, offset)
    },
    re_formula = if (keep_factor_only) ~0 else NULL,
    offset = offset
  )

  data_observed <- fit$data
  data_observed[[factor]] <- as.character(data_observed[[factor]])
  data_observed$is_zero <- data_observed[[abundance]] == 0
  data_observed$normalised <- adjusted[["adjusted___Estimate"]]
  data_observed
}

#' Posterior predictive counts for [plot_boxplot()]
#'
#' Draws from [brms::posterior_predict()], divided by `exp(offset)` so they
#' sit on the same scale as the observed (or adjusted) counts. `factor = NULL`
#' (the default) predicts from the whole model. Naming a factor predicts from
#' that effect alone: group-level effects are dropped by `re_formula`, and the
#' remaining covariates are set to `NA`, which zeroes their contribution to the
#' brms design row. That is the same reduction [adjust_gene()] applies to the
#' counts, so the two stay on the same footing.
#'
#' @param fit A `brmsfit` from [estimate_gene()].
#' @param factor Factor to keep. `NULL` uses all effects.
#' @param offset Name of the offset column.
#' @param number_of_draws Number of posterior predictive draws.
#' @param seed Optional seed passed to [brms::posterior_predict()].
#'
#' @return A tibble of the model data repeated per draw, plus `generated`.
#'
#' @noRd
gene_plot_posterior_predictive_check <- function(fit,
                                                 factor = NULL,
                                                 offset,
                                                 number_of_draws,
                                                 seed = NULL) {
  re_formula <- NULL
  newdata <- fit$data
  if (!is.null(factor)) {
    re_formula <- NA
    newdata <- nullify_newdata(
      fit$data,
      nullify = gene_nuisance_fixed_effects(
        fit,
        factor,
        gene_response_name(fit),
        offset
      ),
      offset = NULL
    )
  }

  predict_args <- list(
    object = fit,
    newdata = newdata,
    ndraws = number_of_draws,
    re_formula = re_formula
  )
  if (!is.null(seed)) {
    predict_args$seed <- seed
  }

  yrep <- do.call(brms::posterior_predict, predict_args)
  yrep <- sweep(yrep, 2, exp(fit$data[[offset]]), FUN = "/")
  n_draw <- nrow(yrep)
  posterior_predictive_check <- tibble::as_tibble(
    fit$data[rep(seq_len(nrow(fit$data)), each = n_draw), , drop = FALSE]
  )
  posterior_predictive_check$generated <- as.vector(yrep)
  posterior_predictive_check
}

#' Boxplot quartiles
#'
#' `fun.data` for the blue posterior predictive boxplots. Whiskers reach the
#' full range of the draws.
#'
#' @param x Numeric vector.
#'
#' @return Named vector of `ymin`, `lower`, `middle`, `upper`, `ymax`.
#'
#' @noRd
calc_boxplot_stat <- function(x) {
  stats <- stats::quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  stats
}

#' Drop a leading zero from axis labels
#'
#' `"0.5"` becomes `".5"`, matching the `sccomp` boxplot scale.
#'
#' @param l Character vector of labels.
#'
#' @return Character vector of the same length.
#'
#' @noRd
drop_leading_zero <- function(l) {
  stringr::str_replace(l, "0(?=.)", "")
}

#' Signed square-root transform
#'
#' `sign(x) * sqrt(abs(x))`, so adjusted counts that go slightly negative
#' still plot.
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector of the same length.
#'
#' @noRd
S_sqrt <- function(x) {
  sign(x) * sqrt(abs(x))
}

#' Inverse of `S_sqrt()`
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector of the same length.
#'
#' @noRd
IS_sqrt <- function(x) {
  x^2 * sign(x)
}

#' ggplot2 transformer for `S_sqrt()`
#'
#' @return A `scales` transform object.
#'
#' @noRd
S_sqrt_trans <- function() {
  scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)
}
