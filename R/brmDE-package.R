#' @keywords internal
"_PACKAGE"

#' @importFrom methods setGeneric setMethod
#' @importFrom rlang .data
#' @importFrom tidytargets tt_iterate tt_data tt_data_list tt_evaluate
#' @importFrom brms hypothesis
NULL

# Symbols resolved by tidytargets::tt_iterate() NSE in the HPC pipeline.
utils::globalVariables(c("brms_fit", "se_input"))
