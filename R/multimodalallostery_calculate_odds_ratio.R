#' calculate odds ratio
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_calculate_odds_ratio()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param case_in_region Input or option used by this workflow; see Usage for the exact interface.
#' @param control_in_region Input or option used by this workflow; see Usage for the exact interface.
#' @param case_out_region Input or option used by this workflow; see Usage for the exact interface.
#' @param control_out_region Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_calculate_odds_ratio <- function (case_in_region, control_in_region, case_out_region, control_out_region) 
{
    contingency_matrix <- matrix(c(case_in_region, control_in_region, case_out_region, control_out_region), nrow = 2, byrow = TRUE)
    fisher_result <- stats::fisher.test(contingency_matrix)
    odds_ratio <- fisher_result$estimate
    p_value <- fisher_result$p.value
    return(list(odds_ratio = odds_ratio, p_value = p_value, matrix = contingency_matrix))
}

