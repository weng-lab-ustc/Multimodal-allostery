#' prepare residue correlations
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_prepare_residue_correlations()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param data Input or option used by this workflow; see Usage for the exact interface.
#' @param alpha Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_prepare_residue_correlations <- function (data, alpha = 0.050000000000000003) 
{
    dplyr::mutate(dplyr::arrange(data, R), residue_order = factor(Pos_real, levels = Pos_real), sig_group = dplyr::case_when(pvalue < 
        alpha & R > 0 ~ "Positive", pvalue < alpha & R < 0 ~ "Negative", TRUE ~ "NS"))
}

