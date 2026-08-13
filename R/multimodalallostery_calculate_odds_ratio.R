#' Calculate Odds Ratio.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param case_in_region Value supplied for `case_in_region`.
#' @param control_in_region Value supplied for `control_in_region`.
#' @param case_out_region Value supplied for `case_out_region`.
#' @param control_out_region Value supplied for `control_out_region`.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_calculate_odds_ratio <- function(case_in_region, control_in_region, case_out_region, control_out_region) {
    contingency_matrix <- matrix(c(case_in_region, control_in_region, case_out_region, control_out_region), nrow = 2, byrow = TRUE)
    fisher_result <- fisher.test(contingency_matrix)
    odds_ratio <- fisher_result$estimate
    p_value <- fisher_result$p.value
    return(list(odds_ratio = odds_ratio, p_value = p_value, matrix = contingency_matrix))
}

