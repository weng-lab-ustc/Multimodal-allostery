#' plot binding interface residue median dd g heatmap
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_binding_interface_residue_median_dd_g_heatmap()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param ddG_file Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_sites Input or option used by this workflow; see Usage for the exact interface.
#' @param position_labels Input or option used by this workflow; see Usage for the exact interface.
#' @param title Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_binding_interface_residue_median_dd_g_heatmap <- function (ddG_file, binding_sites, position_labels, title = "Binding Interface Residues - Median \u0394\u0394Gb") 
{
    ddG <- .ma_read_table(ddG_file)
    ddG <- ddG[, `:=`(Pos_real, Pos + 1)]
    median_values <- ddG %>% dplyr::filter(Pos_real %in% binding_sites) %>% dplyr::group_by(Pos_real) %>% dplyr::summarise(median_ddG = stats::median(`mean_kcal/mol`, 
        na.rm = TRUE)) %>% dplyr::ungroup() %>% dplyr::mutate(residue_label = factor(position_labels[as.character(Pos_real)], 
        levels = position_labels), ) %>% dplyr::arrange(residue_label)
    p <- ggplot2::ggplot(median_values, ggplot2::aes(x = 1, y = residue_label, fill = median_ddG)) + ggplot2::geom_tile(color = "white", 
        linewidth = 0.10000000000000001) + ggplot2::scale_fill_gradient2(low = "#1B38A6", mid = "gray", high = "#F4270C", 
        midpoint = 0, limits = c(-1, 2.5), name = expression(Delta * Delta * "Gb (kcal/mol)")) + ggplot2::labs(title = title, 
        x = NULL, y = "Residue") + ggplot2::theme_minimal() + ggplot2::theme(text = ggplot2::element_text(size = 8), axis.text.x = ggplot2::element_blank(), 
        axis.ticks.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5, 
            size = 8), axis.text.y = ggplot2::element_text(size = 8), axis.title.y = ggplot2::element_text(size = 8), legend.position = "right", 
        panel.border = ggplot2::element_rect(color = "gray90", fill = NA, linewidth = 0.80000000000000004), legend.title = ggplot2::element_text(size = 8), 
        legend.text = ggplot2::element_text(size = 8)) + ggplot2::scale_y_discrete(limits = rev(levels(median_values$residue_label))) + 
        ggplot2::coord_fixed(ratio = 0.40000000000000002)
    return(p)
}

