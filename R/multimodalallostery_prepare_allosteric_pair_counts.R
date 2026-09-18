#' prepare allosteric pair counts
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_prepare_allosteric_pair_counts()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param data Input or option used by this workflow; see Usage for the exact interface.
#' @param pair_order Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_prepare_allosteric_pair_counts <- function (data, pair_order) 
{
    dplyr::mutate(tidyr::pivot_longer(dplyr::select(dplyr::mutate(dplyr::filter(data, Pair %in% pair_order), Pair = factor(Pair, 
        levels = pair_order)), Pair, Allosteric_only_in_X, Allosteric_only_in_Y), cols = c(Allosteric_only_in_X, Allosteric_only_in_Y), 
        names_to = "Type", values_to = "Count"), Type = factor(Type, levels = c("Allosteric_only_in_X", "Allosteric_only_in_Y"), 
        labels = c("Allosteric only in X", "Allosteric only in Y")))
}

