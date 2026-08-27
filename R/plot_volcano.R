#' Volcano plot of a gene-wise pipeline
#'
#' Effect size against null probability for every gene, in the same shape as
#' the `tidybulk` volcano: `log2_fold_change` on x, `pH0` on a reversed log
#' scale on y, so the genes a run is about sit at the top.
#'
#' This takes the pipeline rather than the table it produces, because the
#' pipeline knows both halves of the plot. The hypothesis table is in its
#' store, and the draws its fits kept are on it as metadata, recorded by
#' [estimate()]. That second number is the one a table cannot supply: `pH0` is
#' counted from draws, so it cannot resolve below `1 / ndraws`, and at
#' transcriptome scale many genes come back as exactly 0 — which a log scale
#' has nowhere to put.
#'
#' Those genes are drawn in a shaded band one decade below the resolution and
#' jittered across it. The yellow square in the legend is that bound,
#' `pH0 < 1e-3` (the resolution, in scientific notation): smaller than these
#' draws can measure. The spread
#' within the band carries no information beyond separating the points, so
#' the axis is left unlabelled there. The dashed line is the resolution
#' itself.
#'
#' The result is a `ggplot`, so a pipeline holding several contrasts can be
#' split with `+ facet_wrap(~ hypothesis)` in the same way the `tidybulk`
#' volcano facets by method.
#'
#' @param .data A pipeline from [brmDE()] with [estimate()] and
#'   [brms::hypothesis()] on it. It is evaluated here if it has not been
#'   already, as printing it would evaluate it; targets skips whatever is
#'   built, so plotting a finished run costs nothing.
#' @param effect Name of the effect size column, `"log2_fold_change"` by
#'   default.
#' @param probability Name of the null probability column, `"pH0"` by default.
#'   Pass `"fdr"` to plot the false discovery rate instead.
#' @param significance_threshold Genes below this are coloured and enlarged.
#'   Defaults to `0.05`.
#' @param seed Seed for the jitter, passed to [ggplot2::position_jitter()].
#'   `NA` (default) gives a different jitter each time.
#'
#' @return A `ggplot` object.
#'
#' @seealso [hypothesis_gene()], [plot_boxplot()]
#'
#' @examples
#' \dontrun{
#' data("airway", package = "airway")
#'
#' se <- airway
#' se$offset <- log(colSums(SummarizedExperiment::assay(se, "counts")))
#' se <- tidybulk::estimate_dispersion(se, formula_abundance = ~ dex + cell)
#' se <- estimate_dispersion_log_sd(se, formula_abundance = ~ dex + cell)
#'
#' pipeline <- se |>
#'   brmDE() |>
#'   estimate(
#'     ~ dex + (1 | cell),
#'     offset = "offset",
#'     dispersion_prior_log_mean = "dispersion_prior_log_mean",
#'     dispersion_prior_log_sd = "dispersion_prior_log_sd",
#'     family = brms::negbinomial()
#'   ) |>
#'   hypothesis("dextrt")
#'
#' # Plotting runs the pipeline if printing it has not already.
#' pipeline |> plot_volcano()
#' }
#'
#' @export
plot_volcano <- function(.data,
                         effect = "log2_fold_change",
                         probability = "pH0",
                         significance_threshold = 0.05,
                         seed = NA) {
  if (!inherits(.data, "brmDE_hpc")) {
    stop(
      "plot_volcano() expects a pipeline from brmDE(), which carries both ",
      "the hypothesis table and the draws it was counted from.",
      call. = FALSE
    )
  }
  if (!"hypothesis_tbl" %in% names(.data)) {
    stop(
      "plot_volcano() needs hypothesis() on the pipeline first ",
      "(missing target 'hypothesis_tbl').",
      call. = FALSE
    )
  }
  # A pipeline is a graph until something asks it for values, and this is
  # asking, so it is evaluated here exactly as printing it would evaluate it.
  # targets skips whatever is already built, so plotting a run again is cheap.
  settings <- tidytargets::tt_metadata(.data)
  volcano_ggplot(
    tidytargets::tt_evaluate(.data),
    number_of_draws = settings$chains * settings$draws_sampling,
    effect = effect,
    probability = probability,
    significance_threshold = significance_threshold,
    seed = seed
  )
}

#' Volcano plot of a hypothesis table
#'
#' The plot itself, for a table and the draws behind it. [plot_volcano()] reads
#' both off a pipeline; this is where a table from anywhere else would go in.
#'
#' @param hypotheses A hypothesis table, as `collect_brmde_hpc()` returns it.
#' @param number_of_draws Post-warmup draws behind the fits, which set the
#'   resolution `1 / number_of_draws`. `NULL` infers it from the smallest
#'   non-zero probability in the table, which is the value a counted `pH0`
#'   attains.
#' @param effect Name of the effect size column.
#' @param probability Name of the null probability column.
#' @param significance_threshold Genes below this are coloured and enlarged.
#' @param seed Seed for the jitter.
#'
#' @return A `ggplot` object.
#'
#' @noRd
volcano_ggplot <- function(hypotheses,
                           number_of_draws,
                           effect = "log2_fold_change",
                           probability = "pH0",
                           significance_threshold = 0.05,
                           seed = NA) {
  probabilities <- hypotheses[[probability]]
  resolution <- volcano_resolution(probabilities, number_of_draws)

  volcano_data <- hypotheses
  volcano_data$unresolved___ <-
    !is.na(probabilities) & probabilities < resolution
  volcano_data$significant___ <- probabilities < significance_threshold
  # Genes the draws could not resolve are parked at the centre of the band and
  # spread over it by the jitter, since log10(0) has nowhere to go.
  volcano_data$probability___ <- ifelse(
    volcano_data$unresolved___,
    resolution / sqrt(10),
    probabilities
  )

  significance <- sprintf("%s < %s", probability, significance_threshold)
  unresolved_label <- sprintf("%s < %.0e", probability, resolution)

  ggplot2::ggplot(
    volcano_data,
    ggplot2::aes(
      x = .data[[effect]],
      y = .data$probability___,
      colour = .data$significant___,
      size = .data$significant___
    )
  ) +
    ggplot2::geom_rect(
      data = data.frame(unresolved = TRUE),
      ggplot2::aes(
        xmin = -Inf,
        xmax = Inf,
        ymin = resolution / 10,
        ymax = resolution,
        fill = .data$unresolved
      ),
      colour = NA,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(
      yintercept = resolution,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      data = volcano_data[!volcano_data$unresolved___, , drop = FALSE]
    ) +
    ggplot2::geom_point(
      data = volcano_data[volcano_data$unresolved___, , drop = FALSE],
      # The band is a decade tall and the jitter is in decades, so 0.4 keeps
      # these points off both edges: none of them lands on the resolution line.
      position = ggplot2::position_jitter(height = 0.4, width = 0, seed = seed)
    ) +
    ggplot2::scale_y_continuous(
      transform = scales::transform_compose("log10", "reverse"),
      # Breaks stop at the resolution: a tick inside the band would read as a
      # probability the draws never measured.
      breaks = 10^(0:ceiling(log10(resolution)))
    ) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "red", `FALSE` = "black"),
      name = significance
    ) +
    ggplot2::scale_size_manual(
      values = c(`TRUE` = 1.5, `FALSE` = 0.5),
      name = significance
    ) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "lightyellow"),
      name = unresolved_label
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(colour = "grey40", linewidth = 0.3)
      )
    ) +
    ggplot2::labs(x = effect, y = probability) +
    ggplot2::theme_bw()
}

#' Smallest null probability the draws can resolve
#'
#' `1 / number_of_draws` when the draws are known. Otherwise the smallest
#' non-zero probability in the table, which is the value a counted `pH0`
#' attains when a single draw falls on the null side.
#'
#' @param probability Numeric vector of null probabilities.
#' @param number_of_draws Post-warmup draws, or `NULL` to infer the resolution
#'   from `probability`.
#'
#' @return A single positive number.
#'
#' @noRd
volcano_resolution <- function(probability, number_of_draws) {
  if (!is.null(number_of_draws)) {
    return(1 / number_of_draws)
  }
  resolved <- probability[!is.na(probability) & probability > 0]
  if (length(resolved) == 0L) {
    stop(
      "Every probability is 0, so the draw resolution cannot be read off ",
      "them. Pass `number_of_draws`.",
      call. = FALSE
    )
  }
  min(resolved)
}
