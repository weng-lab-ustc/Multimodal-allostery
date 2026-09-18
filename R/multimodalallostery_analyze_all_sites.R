#' analyze all sites
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_analyze_all_sites()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input_x Input or option used by this workflow; see Usage for the exact interface.
#' @param input_y Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_x Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_y Input or option used by this workflow; see Usage for the exact interface.
#' @param anno Input or option used by this workflow; see Usage for the exact interface.
#' @param fixed_threshold Input or option used by this workflow; see Usage for the exact interface.
#' @param legend_order Input or option used by this workflow; see Usage for the exact interface.
#' @param nbp_residues Input or option used by this workflow; see Usage for the exact interface.
#' @param switch_ii_residues Input or option used by this workflow; see Usage for the exact interface.
#' @param switch_i_residues Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_analyze_all_sites <- function (input_x, input_y, assay_x, assay_y, anno, fixed_threshold = 0.40000000000000002, legend_order = c("Both promoting", 
    "Both disrupting", "Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y", "Allosteric only in X", "Allosteric only in Y", 
    "Other (neutral in both)", "Not significant (FDR >= 0.05)"), nbp_residues = c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 
    32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147), switch_ii_residues = c(58, 59, 60, 61, 62, 63, 64, 65, 66, 
    67, 68, 69, 70, 71, 72, 73, 74, 75, 76), switch_i_residues = c(25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 
    39, 40)) 
{
    prepared <- multimodalallostery_prepare_site_merged_data_with_fdr(input_x = input_x, input_y = input_y, assay_x = assay_x, 
        assay_y = assay_y, anno = anno, fixed_threshold = fixed_threshold)
    merged_data <- prepared$data
    threshold_x <- prepared$threshold_x
    threshold_y <- prepared$threshold_y
    if (!"mt_codon" %in% colnames(merged_data)) {
        merged_data[, `:=`(mt_codon, ifelse(grepl("^WT.*WT$", mt), "WT", substr(mt, nchar(mt), nchar(mt))))]
    }
    merged_data[, `:=`(pass_FDR_x = NA, pass_FDR_y = NA)]
    merged_data[mt_codon == "WT", `:=`(pass_FDR_x = NA, pass_FDR_y = NA, p_adj_x = NA, p_adj_y = NA)]
    merged_data[mt_codon != "WT", `:=`(pass_FDR_x = p_adj_x < 0.050000000000000003, pass_FDR_y = p_adj_y < 0.050000000000000003)]
    merged_data[, `:=`(direction_class, NA_character_)]
    merged_data[mt_codon == "WT", `:=`(direction_class, "WT")]
    merged_data[mt_codon != "WT", `:=`(direction_class, multimodalallostery_classify_by_direction(ddG_x = get(paste0("ddG_", 
        assay_x)), ddG_y = get(paste0("ddG_", assay_y)), threshold_x = fixed_threshold, threshold_y = fixed_threshold))]
    merged_data[, `:=`(final_classification, NA_character_)]
    merged_data[mt_codon == "WT", `:=`(final_classification, "WT")]
    merged_data[mt_codon != "WT", `:=`(final_classification, multimodalallostery_reclassify_by_fdr(direction_class = direction_class, 
        pass_FDR_x = pass_FDR_x, pass_FDR_y = pass_FDR_y))]
    legend_order_with_WT <- c("WT", legend_order)
    merged_data[, `:=`(final_classification, factor(final_classification, levels = legend_order_with_WT))]
    merged_data <- multimodalallostery_add_region_marker(data = merged_data, nbp_residues = nbp_residues, switch_ii_residues = switch_ii_residues, 
        switch_i_residues = switch_i_residues)
    merged_data[mt_codon == "WT", `:=`(plot_group, "WT")]
    threshold_vector <- stats::setNames(c(threshold_x, threshold_y), c(assay_x, assay_y))
    return(list(data = merged_data, thresholds = threshold_vector, assays = c(assay_x, assay_y)))
}

