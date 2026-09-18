#' reclassify by fdr
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_reclassify_by_fdr()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param direction_class Input or option used by this workflow; see Usage for the exact interface.
#' @param pass_FDR_x Input or option used by this workflow; see Usage for the exact interface.
#' @param pass_FDR_y Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_reclassify_by_fdr <- function (direction_class, pass_FDR_x, pass_FDR_y) 
{
    result <- direction_class
    for (i in 1:length(direction_class)) {
        if (direction_class[i] %in% c("Both promoting", "Both disrupting", "Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y")) {
            if (!(pass_FDR_x[i] & pass_FDR_y[i])) {
                if (pass_FDR_x[i] & !pass_FDR_y[i]) {
                  result[i] <- "Allosteric only in X"
                }
                else if (!pass_FDR_x[i] & pass_FDR_y[i]) {
                  result[i] <- "Allosteric only in Y"
                }
                else {
                  result[i] <- "Not significant (FDR >= 0.05)"
                }
            }
        }
        else if (direction_class[i] == "Allosteric only in X") {
            if (!pass_FDR_x[i]) {
                result[i] <- "Not significant (FDR >= 0.05)"
            }
        }
        else if (direction_class[i] == "Allosteric only in Y") {
            if (!pass_FDR_y[i]) {
                result[i] <- "Not significant (FDR >= 0.05)"
            }
        }
    }
    return(result)
}

