#' Reclassify by Fdr.
#'
#' Reusable classification/annotation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param direction_class Value supplied for `direction_class`.
#' @param pass_FDR_x Value supplied for `pass_FDR_x`.
#' @param pass_FDR_y Value supplied for `pass_FDR_y`.
#'
#' @return The input data with classifications or annotations added.
#' @export
multimodalallostery_reclassify_by_fdr <- function(direction_class, pass_FDR_x, pass_FDR_y) {
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

