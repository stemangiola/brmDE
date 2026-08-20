#' Add a tidybulk TMM offset to a SummarizedExperiment
#'
#' Runs [tidybulk::scale_abundance()] on the full object (not gene-wise) and
#' stores `offset = log(1 / multiplier)` in `colData`, matching the
#' immuneBodyMap / HPCell convention. `estimate_gene()` does not compute
#' this offset; this helper (and [brmDE()]) does.
#'
#' @param se A `SummarizedExperiment`.
#' @param abundance Assay name passed to tidybulk (default first assay /
#'   `"counts"`).
#' @param method tidybulk scaling method. Default `"TMM"` (HPCell). Use
#'   `"TMMwsp"` to match the immuneBodyMap pipeline.
#'
#' @return `se` with assays `*_scaled`, and `colData` columns `TMM`,
#'   `multiplier`, and `offset`.
#'
#' @examples
#' \dontrun{
#' data("airway", package = "airway")
#' se <- add_tidybulk_offset(airway[1:150, ])
#' }
#'
#' @export
add_tidybulk_offset <- function(se, abundance = "counts", method = "TMM") {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }
  if (!requireNamespace("tidybulk", quietly = TRUE)) {
    stop("Install tidybulk to calculate the offset.", call. = FALSE)
  }

  se <- tidybulk::identify_abundant(se, abundance = abundance)
  se <- tidybulk::scale_abundance(
    se,
    abundance = abundance,
    method = method
  )
  cd <- SummarizedExperiment::colData(se)
  if (!"multiplier" %in% colnames(cd)) {
    stop(
      "tidybulk::scale_abundance() did not add a `multiplier` column.",
      call. = FALSE
    )
  }
  se$offset <- log(1 / se$multiplier)
  se
}
