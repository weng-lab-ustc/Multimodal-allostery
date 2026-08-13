#' Prepare Region Merged Data with Fdr.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input_x Value supplied for `input_x`.
#' @param input_y Value supplied for `input_y`.
#' @param assay_x Identifier for the first assay.
#' @param assay_y Identifier for the second assay.
#' @param anno Annotation data or annotation file path used by the analysis.
#' @param binding_sites_map Named list mapping assays to binding-site residue positions.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_region_merged_data_with_fdr <- function(input_x, input_y, assay_x, assay_y, anno, binding_sites_map) {
    data_x <- multimodalallostery_load_mutation_data(input_x, assay_x)
    data_y <- multimodalallostery_load_mutation_data(input_y, assay_y)
    threshold_x <- multimodalallostery_calculate_threshold(data_x, assay_x, anno)
    threshold_y <- multimodalallostery_calculate_threshold(data_y, assay_y, anno)
    if (!is.null(binding_sites_map[[assay_x]])) {
        data_x <- data_x[!(Pos_real %in% binding_sites_map[[assay_x]])]
    }
    if (!is.null(binding_sites_map[[assay_y]])) {
        data_y <- data_y[!(Pos_real %in% binding_sites_map[[assay_y]])]
    }
    data_x_clean <- data_x[, .(mt, Pos_real, ddG, ddG_std)]
    setnames(data_x_clean, "ddG", paste0("ddG_", assay_x))
    setnames(data_x_clean, "ddG_std", paste0("std_", assay_x))
    data_y_clean <- data_y[, .(mt, Pos_real, ddG, ddG_std)]
    setnames(data_y_clean, "ddG", paste0("ddG_", assay_y))
    setnames(data_y_clean, "ddG_std", paste0("std_", assay_y))
    merged_data <- merge(data_x_clean, data_y_clean, by = c("mt", "Pos_real"))
    multimodalallostery_pvalue_threshold <- function(av, se, threshold) {
        zscore <- (abs(av) - threshold)/se
        2 * pnorm(abs(zscore), lower.tail = FALSE)
    }
    merged_data[, `:=`(p_x, multimodalallostery_pvalue_threshold(get(paste0("ddG_", assay_x)), get(paste0("std_", assay_x)), threshold_x))]
    merged_data[, `:=`(p_y, multimodalallostery_pvalue_threshold(get(paste0("ddG_", assay_y)), get(paste0("std_", assay_y)), threshold_y))]
    merged_data[, `:=`(p_adj_x, p.adjust(p_x, "BH"))]
    merged_data[, `:=`(p_adj_y, p.adjust(p_y, "BH"))]
    list(data = merged_data, threshold_x = threshold_x, threshold_y = threshold_y)
}

