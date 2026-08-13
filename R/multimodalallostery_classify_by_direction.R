#' Classify by Direction.
#'
#' Reusable classification/annotation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param ddG_x Value supplied for `ddG_x`.
#' @param ddG_y Value supplied for `ddG_y`.
#' @param threshold_x Value supplied for `threshold_x`.
#' @param threshold_y Value supplied for `threshold_y`.
#'
#' @return The input data with classifications or annotations added.
#' @export
multimodalallostery_classify_by_direction <- function(ddG_x, ddG_y, threshold_x, threshold_y) {
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

