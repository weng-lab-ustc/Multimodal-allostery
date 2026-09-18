#' classify by direction
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_classify_by_direction()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param ddG_x Input or option used by this workflow; see Usage for the exact interface.
#' @param ddG_y Input or option used by this workflow; see Usage for the exact interface.
#' @param threshold_x Input or option used by this workflow; see Usage for the exact interface.
#' @param threshold_y Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_classify_by_direction <- function (ddG_x, ddG_y, threshold_x, threshold_y) 
{
    sig_x <- abs(ddG_x) > threshold_x
    sig_y <- abs(ddG_y) > threshold_y
    disrupt_x <- ddG_x > threshold_x
    promote_x <- ddG_x < -threshold_x
    disrupt_y <- ddG_y > threshold_y
    promote_y <- ddG_y < -threshold_y
    result <- rep("neutral", length(ddG_x))
    result[sig_x & sig_y & promote_x & promote_y] <- "Both promoting"
    result[sig_x & sig_y & disrupt_x & disrupt_y] <- "Both disrupting"
    result[sig_x & sig_y & promote_x & disrupt_y] <- "Promoting in X / Disrupting in Y"
    result[sig_x & sig_y & disrupt_x & promote_y] <- "Disrupting in X / Promoting in Y"
    result[sig_x & !sig_y] <- "Allosteric only in X"
    result[!sig_x & sig_y] <- "Allosteric only in Y"
    result[!sig_x & !sig_y] <- "Not significant (FDR >= 0.05)"
    return(result)
}

