#' prepare region merged data with fdr
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_prepare_region_merged_data_with_fdr()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input_x Input or option used by this workflow; see Usage for the exact interface.
#' @param input_y Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_x Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_y Input or option used by this workflow; see Usage for the exact interface.
#' @param anno Input or option used by this workflow; see Usage for the exact interface.
#' @param fixed_threshold Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_sites_map Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_prepare_region_merged_data_with_fdr <- function (input_x, input_y, assay_x, assay_y, anno, fixed_threshold = 0.40000000000000002, binding_sites_map = list(RAF1 = c(21, 
    25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 
    73, 74), K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71), K13 = c(63, 68, 87, 88, 90, 91, 
    92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138), K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 
    97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137))) 
{
    data_x <- multimodalallostery_load_mutation_data(input = input_x, assay_sele = assay_x)
    data_y <- multimodalallostery_load_mutation_data(input = input_y, assay_sele = assay_y)
    threshold_x <- fixed_threshold
    threshold_y <- fixed_threshold
    if (!is.null(binding_sites_map[[assay_x]])) {
        data_x <- data_x[!(Pos_real %in% binding_sites_map[[assay_x]])]
    }
    if (!is.null(binding_sites_map[[assay_y]])) {
        data_y <- data_y[!(Pos_real %in% binding_sites_map[[assay_y]])]
    }
    data_x_clean <- data_x[, .(mt, Pos_real, ddG, ddG_std)]
    data.table::setnames(data_x_clean, "ddG", paste0("ddG_", assay_x))
    data.table::setnames(data_x_clean, "ddG_std", paste0("std_", assay_x))
    data_y_clean <- data_y[, .(mt, Pos_real, ddG, ddG_std)]
    data.table::setnames(data_y_clean, "ddG", paste0("ddG_", assay_y))
    data.table::setnames(data_y_clean, "ddG_std", paste0("std_", assay_y))
    merged_data <- merge(data_x_clean, data_y_clean, by = c("mt", "Pos_real"))
    pvalue_threshold <- function(av, se, threshold) {
        zscore <- (abs(av) - threshold)/se
        2 * stats::pnorm(abs(zscore), lower.tail = FALSE)
    }
    merged_data[, `:=`(p_x, pvalue_threshold(get(paste0("ddG_", assay_x)), get(paste0("std_", assay_x)), threshold_x))]
    merged_data[, `:=`(p_y, pvalue_threshold(get(paste0("ddG_", assay_y)), get(paste0("std_", assay_y)), threshold_y))]
    merged_data[, `:=`(p_adj_x, stats::p.adjust(p_x, "BH"))]
    merged_data[, `:=`(p_adj_y, stats::p.adjust(p_y, "BH"))]
    list(data = merged_data, threshold_x = threshold_x, threshold_y = threshold_y)
}

