#' Calculate Correlation Simple.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input_x Value supplied for `input_x`.
#' @param input_y Value supplied for `input_y`.
#' @param assay_x Identifier for the first assay.
#' @param assay_y Identifier for the second assay.
#' @param binding_sites_map Named assay-to-binding-site residue mapping.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_calculate_correlation_simple <- function(input_x, input_y, assay_x, assay_y, binding_sites_map) {
    data_x <- multimodalallostery_load_dd_g_data(input_x, assay_x)
    data_y <- multimodalallostery_load_dd_g_data(input_y, assay_y)
    if (assay_x %in% names(binding_sites_map)) {
        data_x <- data_x[!(Pos_real %in% binding_sites_map[[assay_x]])]
    }
    if (assay_y %in% names(binding_sites_map)) {
        data_y <- data_y[!(Pos_real %in% binding_sites_map[[assay_y]])]
    }
    merged_data <- merge(data_x[, .(mt, Pos_real, ddG)], data_y[, .(mt, Pos_real, ddG)], by = c("mt", "Pos_real"), suffixes = c(paste0("_", 
        assay_x), paste0("_", assay_y)))
    cor_test <- cor.test(merged_data[[paste0("ddG_", assay_x)]], merged_data[[paste0("ddG_", assay_y)]])
    return(list(r = round(cor_test$estimate, 3), p = cor_test$p.value, n = nrow(merged_data)))
}

