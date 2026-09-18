#' plot region enrichment
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_region_enrichment()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param prepared Input or option used by this workflow; see Usage for the exact interface.
#' @param color_map Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_region_enrichment <- function (prepared, color_map) 
{
    ggplot2::ggplot(prepared$plot_data, ggplot2::aes(x = region, y = frac, fill = category)) + ggplot2::geom_col(position = ggplot2::position_dodge(0.80000000000000004), 
        width = 0.69999999999999996) + 
        #ggplot2::geom_errorbar(ggplot2::aes(ymin = frac - se, ymax = frac + se), position = ggplot2::position_dodge(0.80000000000000004), 
        #width = 0.14999999999999999, linewidth = 0.80000000000000004, color = "#F1DD10", alpha = 0.80000000000000004, na.rm = TRUE) + 
        #ggplot2::geom_errorbar(ggplot2::aes(ymin = frac - se, ymax = frac - se), position = ggplot2::position_dodge(0.80000000000000004), 
        #width = 0.29999999999999999, linewidth = 0.80000000000000004, color = "#F1DD10", alpha = 0.80000000000000004,  na.rm = TRUE) + 
        #ggplot2::geom_errorbar(ggplot2::aes(ymin = frac + se, ymax = frac + se), position = ggplot2::position_dodge(0.80000000000000004), 
        #width = 0.29999999999999999, linewidth = 0.80000000000000004, color = "#F1DD10", alpha = 0.80000000000000004, na.rm = TRUE) + 
        ggplot2::geom_text(data = prepared$odds_data, ggplot2::aes(x = as.numeric(region) + x_offset, y = prepared$fixed_y_position, 
            label = label), size = 3.5, angle = 45, hjust = 0.5, vjust = 0, color = ifelse(prepared$odds_data$category == 
            "Correlated", "#F4AD0C", ifelse(prepared$odds_data$category == "Anti-correlated", "#1B38A6", "grey50"))) + ggplot2::scale_fill_manual(values = color_map) + 
        ggplot2::scale_y_continuous(limits = c(0, prepared$max_y), expand = ggplot2::expansion(mult = c(0, 0.02)), breaks = seq(0, 
            1, 0.20000000000000001)) + ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.5, 
        linewidth = 0.80000000000000004) + ggplot2::labs(y = "Fraction of mutations in region", x = "Structural region") + 
        ggplot2::theme_classic(base_size = 15) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, 
        size = 11), axis.text.y = ggplot2::element_text(size = 15), axis.line = ggplot2::element_line(color = "black", linewidth = 0.5), 
        axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.5), panel.grid = ggplot2::element_blank(), legend.title = ggplot2::element_blank(), 
        legend.position = "bottom", plot.margin = ggplot2::margin(t = 40, r = 10, b = 10, l = 10))
}

