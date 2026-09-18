#' get legend plot
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_get_legend_plot()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param analysis_result Input or option used by this workflow; see Usage for the exact interface.
#' @param base_size Input or option used by this workflow; see Usage for the exact interface.
#' @param color_map Input or option used by this workflow; see Usage for the exact interface.
#' @param legend_order Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_get_legend_plot <- function (analysis_result, base_size = 12, color_map = c(`Not significant (FDR >= 0.05)` = "grey90", `Other (neutral in both)` = "grey90", 
    `Both promoting` = "#FFB0A5", `Both disrupting` = "#F4270C", `Promoting in X / Disrupting in Y` = "#F4AD0C", `Disrupting in X / Promoting in Y` = "#F1DD10", 
    `Allosteric only in X` = "#1B38A6", `Allosteric only in Y` = "#75C2F6"), legend_order = c("Both promoting", "Both disrupting", 
    "Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y", "Allosteric only in X", "Allosteric only in Y", 
    "Other (neutral in both)", "Not significant (FDR >= 0.05)")) 
{
    merged_data <- analysis_result$data
    protein_x <- analysis_result$names[1]
    protein_y <- analysis_result$names[2]
    threshold_x <- analysis_result$thresholds[1]
    threshold_y <- analysis_result$thresholds[2]
    ddG_x_col <- paste0("ddG_", protein_x)
    ddG_y_col <- paste0("ddG_", protein_y)
    dummy_data <- data.frame(x = 0, y = 0, final_classification = factor(legend_order, levels = legend_order))
    colnames(dummy_data)[1:2] <- c(ddG_x_col, ddG_y_col)
    p_legend <- ggplot2::ggplot(dummy_data) + ggplot2::geom_point(.ma_aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"), 
        size = 3, stroke = 0.29999999999999999) + ggplot2::scale_color_manual(values = color_map, breaks = legend_order, 
        drop = FALSE) + ggplot2::theme_classic(base_size = base_size) + ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = base_size - 
        1), legend.title = ggplot2::element_blank(), legend.key.size = grid::unit(0.40000000000000002, "cm"), legend.spacing.y = grid::unit(0.050000000000000003, 
        "cm"), legend.margin = ggplot2::margin(t = 5, b = 5), panel.border = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), 
        axis.title = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) + ggplot2::guides(color = ggplot2::guide_legend(ncol = 2, 
        byrow = TRUE, override.aes = list(size = 3)))
    return(p_legend)
}

