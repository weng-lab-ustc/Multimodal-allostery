#' define chemotype
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_define_chemotype()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param aa Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_define_chemotype <- function (aa) 
{
    aromatic <- c("F", "W", "Y")
    aliphatic <- c("A", "V", "I", "L", "M")
    polar_uncharged <- c("S", "T", "N", "Q", "C")
    positive <- c("K", "R", "H")
    negative <- c("D", "E")
    special <- c("G", "P")
    if (aa %in% aromatic) 
        return("Aromatic")
    if (aa %in% aliphatic) 
        return("Aliphatic")
    if (aa %in% polar_uncharged) 
        return("Polar uncharged")
    if (aa %in% positive) 
        return("Positive")
    if (aa %in% negative) 
        return("Negative")
    if (aa %in% special) 
        return("Special")
    return(NA)
}

