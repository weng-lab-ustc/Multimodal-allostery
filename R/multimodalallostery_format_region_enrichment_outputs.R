#' Format Region Enrichment Outputs.
#'
#' Reusable utility logic consolidated from the selected Figure/Panel implementations.
#'
#' @param result Analysis result to summarize or plot.
#'
#' @return The transformed result produced by the function.
#' @export
multimodalallostery_format_region_enrichment_outputs <- function(result) {
    odds <- result$or[, .(pair, region, category, OR, OR_low, OR_high, p)]
    odds[, `:=`(p_signif, ifelse(p < 0.05, "*", "ns"))]
    odds[, `:=`(OR_text, paste0(round(OR, 2), ifelse(p < 0.05, "*", "")))]
    fractions <- result$plot[, .(pair, region, category, frac, se, n)]
    list(odds = odds, fractions = fractions)
}

