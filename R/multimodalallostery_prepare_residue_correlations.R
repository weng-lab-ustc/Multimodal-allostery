#' Prepare Residue Correlations.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data Input data frame or data table.
#' @param alpha Value supplied for `alpha`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_residue_correlations <- function(data, alpha = 0.05) {
    dplyr::mutate(dplyr::arrange(data, R), residue_order = factor(Pos_real, levels = Pos_real), sig_group = dplyr::case_when(pvalue < 
        alpha & R > 0 ~ "Positive", pvalue < alpha & R < 0 ~ "Negative", TRUE ~ "NS"))
}

