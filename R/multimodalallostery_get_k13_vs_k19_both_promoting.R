#' get k13 vs k19 both promoting
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_get_k13_vs_k19_both_promoting()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input_k13 Input or option used by this workflow; see Usage for the exact interface.
#' @param input_k19 Input or option used by this workflow; see Usage for the exact interface.
#' @param anno Input or option used by this workflow; see Usage for the exact interface.
#' @param fixed_threshold Input or option used by this workflow; see Usage for the exact interface.
#' @param nbp_residues Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_sites_map Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_get_k13_vs_k19_both_promoting <- function (input_k13, input_k19, anno, fixed_threshold = 0.40000000000000002, nbp_residues = c(12, 13, 14, 15, 16, 17, 18, 
    28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147), binding_sites_map = list(RAF1 = c(21, 25, 29, 
    31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71), 
    K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138), K19 = c(68, 
        87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137))) 
{
    input_files <- list(K13 = input_k13, K19 = input_k19)
    anno <- anno
    prepared <- multimodalallostery_prepare_mapped_merged_data_with_fdr(input_x = input_files$K13, input_y = input_files$K19, 
        assay_x = "K13", assay_y = "K19", anno = anno, binding_sites_map = binding_sites_map, fixed_threshold = fixed_threshold)
    merged_data <- prepared$data
    threshold_K13 <- prepared$threshold_x
    threshold_K19 <- prepared$threshold_y
    cat("\nFixed threshold used:", fixed_threshold, "kcal/mol\n")
    cat("Threshold K13:", threshold_K13, "\n")
    cat("Threshold K19:", threshold_K19, "\n")
    merged_data[, `:=`(pass_FDR_K13, p_adj_x < 0.050000000000000003)]
    merged_data[, `:=`(pass_FDR_K19, p_adj_y < 0.050000000000000003)]
    merged_data[, `:=`(direction_class, multimodalallostery_classify_by_direction(ddG_x = ddG_K13, ddG_y = ddG_K19, threshold_x = threshold_K13, 
        threshold_y = threshold_K19))]
    merged_data[, `:=`(final_classification, multimodalallostery_reclassify_by_fdr(direction_class = direction_class, pass_FDR_x = pass_FDR_K13, 
        pass_FDR_y = pass_FDR_K19))]
    both_promoting_muts <- merged_data[final_classification == "Both promoting"]
    both_promoting_muts[, `:=`(is_NBP, Pos_real %in% nbp_residues)]
    cat("\nK13 vs K19 Both promoting mutations found:", nrow(both_promoting_muts), "\n")
    cat("  - NBP residues:", sum(both_promoting_muts$is_NBP), "\n")
    cat("  - Non-NBP residues:", sum(!both_promoting_muts$is_NBP), "\n")
    print(both_promoting_muts[, .(mt, Pos_real, is_NBP, ddG_K13, ddG_K19)])
    return(both_promoting_muts)
}

