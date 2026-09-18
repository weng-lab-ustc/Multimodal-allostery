#' calculate correlation simple
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_calculate_correlation_simple()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input_x Input or option used by this workflow; see Usage for the exact interface.
#' @param input_y Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_x Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_y Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_sites_map Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_calculate_correlation_simple <- function (input_x, input_y, assay_x, assay_y, binding_sites_map = list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 
    67, 71), RALGDS = c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67), PI3KCG = c(3, 21, 24, 25, 33, 36, 37, 38, 39, 
    40, 41, 63, 64, 70, 73), SOS1 = c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 
    65, 66, 67, 70, 149, 153), K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74), K27 = c(21, 
    24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71), K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 
    98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138), K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 
    105, 107, 108, 125, 129, 133, 136, 137))) 
{
    data_x <- multimodalallostery_load_dd_g_data(input = input_x, assay_sele = assay_x)
    data_y <- multimodalallostery_load_dd_g_data(input = input_y, assay_sele = assay_y)
    assay_map <- c(RAF1 = "RAF1", RAL = "RALGDS", PI3 = "PI3KCG", SOS = "SOS1", K55 = "K55", K27 = "K27", K13 = "K13", K19 = "K19" )
    site_x <- binding_sites_map[[assay_map[[assay_x]]]]
    site_y <- binding_sites_map[[assay_map[[assay_y]]]]
    data_x <- data_x[!(Pos_real %in% site_x)]
    data_y <- data_y[!(Pos_real %in% site_y)]
    merged_data <- merge(data_x[, .(mt, Pos_real, ddG)], data_y[, .(mt, Pos_real, ddG)], by = c("mt", "Pos_real"), suffixes = c(paste0("_", assay_x), paste0("_", assay_y)))
    cor_test <- stats::cor.test(merged_data[[paste0("ddG_", assay_x)]],merged_data[[paste0("ddG_", assay_y)]])
    list(r = round(unname(cor_test$estimate), 3), p = cor_test$p.value, n = nrow(merged_data))
}

