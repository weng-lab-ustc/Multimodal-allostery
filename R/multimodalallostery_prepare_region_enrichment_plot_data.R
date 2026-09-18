#' prepare region enrichment plot data
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_prepare_region_enrichment_plot_data()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param result Input or option used by this workflow; see Usage for the exact interface.
#' @param target_pair Input or option used by this workflow; see Usage for the exact interface.
#' @param region_levels Input or option used by this workflow; see Usage for the exact interface.
#' @param categories Input or option used by this workflow; see Usage for the exact interface.
#' @param fixed_y_position Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_prepare_region_enrichment_plot_data <- function (result, target_pair, region_levels, categories = c("Correlated", "Anti-correlated", "Other"), fixed_y_position = 1.2) 
{
    plot_data <- data.table::copy(result$plot[category %in% categories])
    plot_data[, `:=`(region, factor(region, levels = region_levels))]
    plot_data[, `:=`(category, factor(category, levels = categories))]
    odds_data <- data.table::copy(result$or[pair == target_pair & category %in% categories])
    odds_data[, `:=`(region, factor(region, levels = region_levels))]
    odds_data[, `:=`(x_offset, dplyr::case_when(category == "Correlated" ~ -0.25, category == "Anti-correlated" ~ 0, TRUE ~ 
        0.25))]
    ##odds_data[, `:=`(label, paste0("OR = ", round(OR, 2), ifelse(p < 0.050000000000000003, "*", "")))]
    odds_data[, `:=`(label, paste0("OR = ", round(OR, 2)))]
    max_y <- max(plot_data$frac + plot_data$se, na.rm = TRUE)
    max_y <- if (fixed_y_position > max_y) 
        fixed_y_position + 0.10000000000000001
    else max_y * 1.1499999999999999
    list(plot_data = plot_data, odds_data = odds_data, fixed_y_position = fixed_y_position, max_y = max(max_y, 1))
}

