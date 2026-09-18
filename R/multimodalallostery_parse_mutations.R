#' parse mutations
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_parse_mutations()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param mut_str Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_parse_mutations <- function (mut_str) 
{
    muts <- trimws(unlist(strsplit(mut_str, ",")))
    unique(data.table::data.table(mutation = unique(muts[muts != ""])))
}

