#' merge dimsum data
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_merge_dimsum_data()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param merge_1 Input or option used by this workflow; see Usage for the exact interface.
#' @param merge_2 Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_merge_dimsum_data <- function (merge_1, merge_2) 
{
    a1 <- as.character(substitute(merge_1))
    a2 <- as.character(substitute(merge_2))
    merge_1[, `:=`(assay, a1)]
    merge_2[, `:=`(assay, a2)]
    output <- rbind(merge_1, merge_2)
    return(output)
}

