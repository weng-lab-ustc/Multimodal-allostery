#' Plot Allosteric Region Enrichment.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param prepared Value supplied for `prepared`.
#' @param color_map Value supplied for `color_map`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_allosteric_region_enrichment <- function(prepared, color_map) {
    ggplot2::ggplot(prepared$plot_data, ggplot2::aes(x = region, y = frac, fill = category)) + ggplot2::geom_col(position = ggplot2::position_dodge(0.9), 
        width = 0.7) + ggplot2::geom_errorbar(ggplot2::aes(ymin = frac - se, ymax = frac + se), position = ggplot2::position_dodge(0.9), 
        width = 0.15, linewidth = 0.8, color = "#F1DD10", alpha = 0.8, na.rm = TRUE) + ggplot2::geom_text(data = prepared$odds_data, 
        ggplot2::aes(x = as.numeric(region) + x_offset, y = prepared$fixed_y_position, label = label, color = category), 
        size = 3.5, angle = 45, hjust = 0.5, vjust = 0) + ggplot2::scale_fill_manual(values = color_map, labels = c(paste0("Allosteric only in ", 
        prepared$assays[1]), paste0("Allosteric only in ", prepared$assays[2]), "Other")) + ggplot2::scale_color_manual(values = color_map, 
        guide = "none") + ggplot2::scale_y_continuous(limits = c(0, prepared$max_y), expand = ggplot2::expansion(mult = c(0, 
        0.02)), breaks = seq(0, 1, 0.2)) + ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.5, 
        linewidth = 0.8) + ggplot2::labs(y = "Fraction of mutations in region", x = "Structural region") + ggplot2::theme_classic(base_size = 15) + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 11), axis.text.y = ggplot2::element_text(size = 15), 
            legend.title = ggplot2::element_blank(), legend.position = "bottom", panel.grid = ggplot2::element_blank(), plot.margin = ggplot2::margin(t = 40, 
                r = 10, b = 10, l = 10))
}

