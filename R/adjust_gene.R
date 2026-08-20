#' Adjust one gene by removing unwanted effects
#'
#' Posterior residuals (optionally divided by `exp(offset)`) are added to
#' fitted values from a newdata grid where selected covariates are set to
#' `NA` and the offset is zeroed. Random effects are included or dropped
#' via `re_formula`. This is the gene-wise form of the immuneBodyMap
#' `remove_unwanted_effect()` step.
#'
#' @param fit A `brmsfit` from [estimate_gene()].
#' @param nullify Character vector of covariate names to set to `NA` in
#'   `newdata` (nuisance terms). Ignored if `newdata` is supplied.
#' @param newdata Optional data frame passed to `brms::fitted()`. If `NULL`,
#'   `fit$data` is used after applying `nullify` and `offset_value`.
#' @param re_formula Random-effect formula for the fitted component.
#'   Use `NA` to drop all random effects; `~(1 | tissue_groups)` keeps
#'   tissue intercepts.
#' @param offset Name of the offset column in `fit$data`.
#' @param offset_value Value written into the offset column of `newdata`
#'   (default `0`).
#' @param robust If `TRUE`, summarise the posterior with medians.
#' @param correct_by_offset If `TRUE`, divide residuals by `exp(offset)`
#'   from the original data.
#' @param sample_id Optional character vector of sample ids to bind onto the
#'   result. A warning is issued if the length does not match.
#'
#' @return A tibble with `adjusted___`, `residuals___`, and `fitted___`
#'   posterior summaries, one row per observation.
#'
#' @examples
#' \dontrun{
#' data("airway", package = "airway")
#' se <- airway["ENSG00000120129", ]
#' se$offset <- log(colSums(SummarizedExperiment::assay(airway, "counts")))
#' se$dex <- relevel(factor(se$dex), ref = "untrt")
#' fit <- estimate_gene(se, ~ dex + (1 | cell), offset = "offset")
#'
#' adjust_gene(
#'   fit,
#'   nullify = "dex",
#'   re_formula = ~(1 | cell)
#' )
#' }
#'
#' @export
adjust_gene <- function(fit,
                        nullify = NULL,
                        newdata = NULL,
                        re_formula = ~0,
                        offset = "offset",
                        offset_value = 0,
                        robust = TRUE,
                        correct_by_offset = TRUE,
                        sample_id = NULL) {
  if (is.null(newdata)) {
    newdata <- nullify_newdata(
      fit$data,
      nullify = nullify,
      offset = offset,
      offset_value = offset_value
    )
  }

  fitted_residuals <- stats::residuals(fit, robust = robust, summary = FALSE)

  if (isTRUE(correct_by_offset)) {
    if (!offset %in% names(fit$data)) {
      stop("Offset column '", offset, "' was not found in `fit$data`.", call. = FALSE)
    }
    fitted_residuals <- sweep(fitted_residuals, 2, exp(fit$data[[offset]]), FUN = "/")
  }

  fitted_values <- stats::fitted(
    fit,
    newdata = newdata,
    re_formula = re_formula,
    summary = FALSE,
    offset = 0
  )

  if (!identical(dim(fitted_residuals), dim(fitted_values))) {
    stop(
      "Residual and fitted draws have different dimensions (",
      paste(dim(fitted_residuals), collapse = " x "), " vs ",
      paste(dim(fitted_values), collapse = " x "),
      "). Adding them would recycle draws across samples.",
      call. = FALSE
    )
  }

  adjusted_counts <- fitted_values + fitted_residuals

  summarise_draws <- function(x, prefix) {
    tbl <- tibble::as_tibble(brms::posterior_summary(x, robust = robust))
    names(tbl) <- paste0(prefix, names(tbl))
    tbl
  }

  out <- dplyr::bind_cols(
    summarise_draws(adjusted_counts, "adjusted___"),
    summarise_draws(fitted_residuals, "residuals___"),
    summarise_draws(fitted_values, "fitted___")
  )

  if (!is.null(sample_id)) {
    if (length(sample_id) != nrow(out)) {
      stop(
        "`sample_id` length (", length(sample_id),
        ") does not match the number of adjusted rows (", nrow(out), ").",
        call. = FALSE
      )
    }
    out$sample_id <- sample_id
  }

  out
}
