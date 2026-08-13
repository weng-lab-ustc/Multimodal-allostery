#' Analyze all Sites.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input_x Value supplied for `input_x`.
#' @param input_y Value supplied for `input_y`.
#' @param assay_x Identifier for the first assay.
#' @param assay_y Identifier for the second assay.
#' @param anno Annotation data or annotation file path used by the analysis.
#' @param legend_order Ordered classification labels used in downstream plots.
#' @param nbp_residues Residue positions assigned to the nucleotide-binding pocket.
#' @param switch_i_residues Residue positions assigned to Switch I.
#' @param switch_ii_residues Residue positions assigned to Switch II.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_analyze_all_sites <- function(input_x, input_y, assay_x, assay_y, anno, legend_order,
                              nbp_residues, switch_i_residues, switch_ii_residues) {
    prepared <- multimodalallostery_prepare_site_merged_data_with_fdr(input_x, input_y, assay_x, assay_y, anno)
    merged_data <- prepared$data
    threshold_x <- prepared$threshold_x
    threshold_y <- prepared$threshold_y
    if (!"mt_codon" %in% colnames(merged_data)) {
        merged_data[, `:=`(mt_codon, ifelse(grepl("^WT.*WT$", mt), "WT", substr(mt, nchar(mt), nchar(mt))))]
    }
    merged_data[, `:=`(pass_FDR_x = NA, pass_FDR_y = NA)]
    merged_data[mt_codon == "WT", `:=`(pass_FDR_x = NA, pass_FDR_y = NA, p_adj_x = NA, p_adj_y = NA)]
    merged_data[mt_codon != "WT", `:=`(pass_FDR_x = p_adj_x < 0.05, pass_FDR_y = p_adj_y < 0.05)]
    merged_data[, `:=`(direction_class, NA_character_)]
    merged_data[mt_codon == "WT", `:=`(direction_class, "WT")]
    merged_data[mt_codon != "WT", `:=`(direction_class, multimodalallostery_classify_by_direction(get(paste0("ddG_", assay_x)), get(paste0("ddG_", 
        assay_y)), threshold_x, threshold_y))]
    merged_data[, `:=`(final_classification, NA_character_)]
    merged_data[mt_codon == "WT", `:=`(final_classification, "WT")]
    merged_data[mt_codon != "WT", `:=`(final_classification, multimodalallostery_reclassify_by_fdr(direction_class, pass_FDR_x, pass_FDR_y))]
    legend_order_with_WT <- c("WT", legend_order)
    merged_data[, `:=`(final_classification, factor(final_classification, levels = legend_order_with_WT))]
    merged_data <- multimodalallostery_add_region_marker(merged_data, nbp_residues, switch_i_residues, switch_ii_residues)
    merged_data[mt_codon == "WT", `:=`(plot_group, "WT")]
    threshold_vector <- setNames(c(threshold_x, threshold_y), c(assay_x, assay_y))
    return(list(data = merged_data, thresholds = threshold_vector, assays = c(assay_x, assay_y)))
}

