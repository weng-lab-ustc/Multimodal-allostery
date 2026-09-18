#' classify two step
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_classify_two_step()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param merged_data Input or option used by this workflow; see Usage for the exact interface.
#' @param threshold_x Input or option used by this workflow; see Usage for the exact interface.
#' @param threshold_y Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_x Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_y Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_classify_two_step <- function (merged_data, threshold_x, threshold_y, assay_x, assay_y) 
{
    ddG_x_col <- paste0("ddG_", assay_x)
    ddG_y_col <- paste0("ddG_", assay_y)
    merged_data[, `:=`(direction_class, multimodalallostery_classify_by_direction(ddG_x = get(ddG_x_col), ddG_y = get(ddG_y_col), 
        threshold_x = threshold_x, threshold_y = threshold_y))]
    merged_data[, `:=`(pass_FDR_x, p_adj_x < 0.050000000000000003)]
    merged_data[, `:=`(pass_FDR_y, p_adj_y < 0.050000000000000003)]
    merged_data[, `:=`(final_classification, multimodalallostery_reclassify_by_fdr(direction_class = direction_class, pass_FDR_x = pass_FDR_x, 
        pass_FDR_y = pass_FDR_y))]
    return(merged_data)
}

