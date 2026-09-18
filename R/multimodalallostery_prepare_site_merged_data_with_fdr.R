#' prepare site merged data with fdr
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_prepare_site_merged_data_with_fdr()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input_x Input or option used by this workflow; see Usage for the exact interface.
#' @param input_y Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_x Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_y Input or option used by this workflow; see Usage for the exact interface.
#' @param anno Input or option used by this workflow; see Usage for the exact interface.
#' @param fixed_threshold Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_prepare_site_merged_data_with_fdr <- function (input_x, input_y, assay_x, assay_y, anno, fixed_threshold = 0.40000000000000002) 
{
    data_x <- multimodalallostery_load_site_mutation_data(input = input_x, assay_sele = assay_x)
    data_y <- multimodalallostery_load_site_mutation_data(input = input_y, assay_sele = assay_y)
    threshold_x <- fixed_threshold
    threshold_y <- fixed_threshold
    data_x_clean <- data_x[, .(mt, Pos_real, mt_codon, ddG, ddG_std)]
    data.table::setnames(data_x_clean, "ddG", paste0("ddG_", assay_x))
    data.table::setnames(data_x_clean, "ddG_std", paste0("std_", assay_x))
    data_y_clean <- data_y[, .(mt, Pos_real, mt_codon, ddG, ddG_std)]
    data.table::setnames(data_y_clean, "ddG", paste0("ddG_", assay_y))
    data.table::setnames(data_y_clean, "ddG_std", paste0("std_", assay_y))
    merged_data <- merge(data_x_clean, data_y_clean, by = c("mt", "Pos_real", "mt_codon"), all = TRUE)
    merged_data[is.na(get(paste0("ddG_", assay_x))), `:=`((paste0("ddG_", assay_x)), 0)]
    merged_data[is.na(get(paste0("ddG_", assay_y))), `:=`((paste0("ddG_", assay_y)), 0)]
    merged_data[is.na(get(paste0("std_", assay_x))), `:=`((paste0("std_", assay_x)), 0)]
    merged_data[is.na(get(paste0("std_", assay_y))), `:=`((paste0("std_", assay_y)), 0)]
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

