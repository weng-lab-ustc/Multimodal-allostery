#' Prepare Allosteric Region Plot Data.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param result Analysis result to summarize or plot.
#' @param target_pair Value supplied for `target_pair`.
#' @param base_regions Value supplied for `base_regions`.
#' @param categories Value supplied for `categories`.
#' @param fixed_y_position Value supplied for `fixed_y_position`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_allosteric_region_plot_data <- function(result, target_pair, base_regions, categories = c("Allosteric only in X", 
    "Allosteric only in Y", "Other"), fixed_y_position = 1.2) {
    plot_data <- data.table::copy(result$plot[category %in% categories])
    region_levels <- c(base_regions, sort(grep("Second Shell", unique(plot_data$region), value = TRUE)))
    plot_data[, `:=`(region, factor(region, levels = region_levels))]
    plot_data[, `:=`(category, factor(category, levels = categories))]
    odds_data <- data.table::copy(result$or[pair == target_pair & category %in% categories & region %in% region_levels])
    odds_data[, `:=`(region, factor(region, levels = region_levels))]
    odds_data[, `:=`(x_offset, dplyr::case_when(category == "Allosteric only in X" ~ -0.25, category == "Allosteric only in Y" ~ 
        0, TRUE ~ 0.25))]
    odds_data[, `:=`(label, paste0("OR = ", round(OR, 2), ifelse(p < 0.05, "*", "")))]
    max_y <- max(plot_data$frac + plot_data$se, na.rm = TRUE)
    max_y <- if (fixed_y_position > max_y) 
        fixed_y_position + 0.1
    else max_y * 1.15
    list(plot_data = plot_data, odds_data = odds_data, region_levels = region_levels, fixed_y_position = fixed_y_position, 
        max_y = max(max_y, 1), assays = strsplit(target_pair, " vs ")[[1]])
}

