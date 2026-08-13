#' Format Allosteric Region Outputs.
#'
#' Reusable utility logic consolidated from the selected Figure/Panel implementations.
#'
#' @param result Analysis result to summarize or plot.
#'
#' @return The transformed result produced by the function.
#' @export
multimodalallostery_format_allosteric_region_outputs <- function(result) {
    output <- result$or[, .(pair, region, category, OR, OR_low, OR_high, p)]
    output[, `:=`(p_signif, ifelse(p < 0.05, "*", "ns"))]
    output[, `:=`(OR_text, paste0(round(OR, 2), ifelse(p < 0.05, "*", "")))]
    output
}

