#' Classify Two Step.
#'
#' Reusable classification/annotation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param merged_data Value supplied for `merged_data`.
#' @param threshold_x Value supplied for `threshold_x`.
#' @param threshold_y Value supplied for `threshold_y`.
#' @param assay_x Identifier for the first assay.
#' @param assay_y Identifier for the second assay.
#'
#' @return The input data with classifications or annotations added.
#' @export
multimodalallostery_classify_two_step <- function(merged_data, threshold_x, threshold_y, assay_x, assay_y) {
    ddG_x_col <- paste0("ddG_", assay_x)
    ddG_y_col <- paste0("ddG_", assay_y)
    merged_data[, `:=`(direction_class, multimodalallostery_classify_by_direction(get(ddG_x_col), get(ddG_y_col), threshold_x, threshold_y))]
    merged_data[, `:=`(pass_FDR_x, p_adj_x < 0.05)]
    merged_data[, `:=`(pass_FDR_y, p_adj_y < 0.05)]
    merged_data[, `:=`(final_classification, multimodalallostery_reclassify_by_fdr(direction_class, pass_FDR_x, pass_FDR_y))]
    return(merged_data)
}

