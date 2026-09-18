#' load k27 vs k19 data
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_load_k27_vs_k19_data()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input_k27 Input or option used by this workflow; see Usage for the exact interface.
#' @param input_k19 Input or option used by this workflow; see Usage for the exact interface.
#' @param anno Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_sites_map Input or option used by this workflow; see Usage for the exact interface.
#' @param fixed_threshold Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_load_k27_vs_k19_data <- function (input_k27, input_k19, anno, binding_sites_map = list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), 
    K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71), K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 
        95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138), K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 
        98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)), fixed_threshold = 0.40000000000000002) 
{
    input_files <- list(K27 = input_k27, K19 = input_k19)
    anno <- anno
    prepared <- multimodalallostery_prepare_mapped_merged_data_with_fdr(input_x = input_files$K27, input_y = input_files$K19, 
        assay_x = "K27", assay_y = "K19", anno = anno, binding_sites_map = binding_sites_map, fixed_threshold = fixed_threshold)
    cat("\nK27 vs K19 thresholds:\n")
    cat("  K27 threshold:", prepared$threshold_x, "\n")
    cat("  K19 threshold:", prepared$threshold_y, "\n")
    return(list(data = prepared$data, threshold_K27 = prepared$threshold_x, threshold_K19 = prepared$threshold_y))
}

