#' Prepare Region Enrichment Plot Data.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param result Analysis result to summarize or plot.
#' @param target_pair Value supplied for `target_pair`.
#' @param region_levels Value supplied for `region_levels`.
#' @param categories Value supplied for `categories`.
#' @param fixed_y_position Value supplied for `fixed_y_position`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_region_enrichment_plot_data <- function(result, target_pair, region_levels, categories = c("Correlated", "Anti-correlated", 
    "Other"), fixed_y_position = 1.2) {
    plot_data <- data.table::copy(result$plot[category %in% categories])
    plot_data[, `:=`(region, factor(region, levels = region_levels))]
    plot_data[, `:=`(category, factor(category, levels = categories))]
    odds_data <- data.table::copy(result$or[pair == target_pair & category %in% categories])
    odds_data[, `:=`(region, factor(region, levels = region_levels))]
    odds_data[, `:=`(x_offset, dplyr::case_when(category == "Correlated" ~ -0.25, category == "Anti-correlated" ~ 0, TRUE ~ 
        0.25))]
    odds_data[, `:=`(label, paste0("OR = ", round(OR, 2)))]
    max_y <- max(plot_data$frac + plot_data$se, na.rm = TRUE)
    max_y <- if (fixed_y_position > max_y) 
        fixed_y_position + 0.1
    else max_y * 1.15
    list(plot_data = plot_data, odds_data = odds_data, fixed_y_position = fixed_y_position, max_y = max(max_y, 1))
}

