#' prepare mutation sets
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_prepare_mutation_sets()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param named_tables Input or option used by this workflow; see Usage for the exact interface.
#' @param mutation_column Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_prepare_mutation_sets <- function (named_tables, mutation_column = "mutation") 
{
    lapply(named_tables, function(x) x[[mutation_column]])
}

