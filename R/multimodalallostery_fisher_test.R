#' fisher test
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_fisher_test()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param case_in Input or option used by this workflow; see Usage for the exact interface.
#' @param case_out Input or option used by this workflow; see Usage for the exact interface.
#' @param control_in Input or option used by this workflow; see Usage for the exact interface.
#' @param control_out Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_fisher_test <- function (case_in, case_out, control_in, control_out) 
{
    mat <- matrix(c(case_in, control_in, case_out, control_out), nrow = 2)
    res <- stats::fisher.test(mat)
    return(list(OR = res$estimate, p = res$p.value))
}

