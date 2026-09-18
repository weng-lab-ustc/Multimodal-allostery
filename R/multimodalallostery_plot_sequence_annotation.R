#' plot sequence annotation
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_sequence_annotation()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param prepared_sequence Input or option used by this workflow; see Usage for the exact interface.
#' @param gtp_pocket Input or option used by this workflow; see Usage for the exact interface.
#' @param functional_loops Input or option used by this workflow; see Usage for the exact interface.
#' @param core_residues Input or option used by this workflow; see Usage for the exact interface.
#' @param beta_sheets Input or option used by this workflow; see Usage for the exact interface.
#' @param alpha_helices Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_sequence_annotation <- function (prepared_sequence, gtp_pocket, functional_loops, core_residues, beta_sheets, alpha_helices) 
{
    sequence_data <- prepared_sequence$sequence_data
    sequence_length <- prepared_sequence$sequence_length
    points_per_mm <- 72.269999999999996/25.399999999999999
    p <- ggplot2::ggplot() + ggplot2::geom_text(data = sequence_data, ggplot2::aes(x = position, y = 0, label = residue, 
        color = binding_group), size = 8/points_per_mm, family = "mono", angle = 90, vjust = 0.5, hjust = 0.5) + ggplot2::scale_color_manual(values = prepared_sequence$binding_colors, 
        labels = c("Other", "K13 binding", "RAF1 binding", "Both"), name = "Binding Sites") + ggplot2::geom_rect(data = data.frame(position = gtp_pocket), 
        ggplot2::aes(xmin = position - 0.5, xmax = position + 0.5, ymin = -0.29999999999999999, ymax = -0.20000000000000001), 
        fill = "#F1DD10", alpha = 0.5) + ggplot2::annotate("text", x = mean(range(gtp_pocket)), y = -0.25, label = "GTP Pocket", 
        color = "#F1DD10", size = 8/points_per_mm) + ggplot2::geom_rect(data = functional_loops, ggplot2::aes(xmin = xstart - 
        0.5, xmax = xend + 0.5, ymin = -0.5, ymax = -0.40000000000000002, fill = col), alpha = 0.29999999999999999) + ggplot2::geom_text(data = functional_loops, 
        ggplot2::aes(x = (xstart + xend)/2, y = -0.45000000000000001, label = col), size = 8/points_per_mm, color = "black") + 
        ggplot2::scale_fill_manual(values = c("#FFB0A5", "#FF0066", "#007A20"), name = "Functional Loops") + ggplot2::geom_rect(data = data.frame(position = core_residues), 
        ggplot2::aes(xmin = position - 0.5, xmax = position + 0.5, ymin = -0.69999999999999996, ymax = -0.59999999999999998), 
        fill = "#85409D", alpha = 0.90000000000000002) + ggplot2::annotate("text", x = mean(range(core_residues)), y = -0.65000000000000002, 
        label = "Core", color = "black", size = 8/points_per_mm) + ggplot2::geom_rect(data = beta_sheets, ggplot2::aes(xmin = xstart - 
        0.5, xmax = xend + 0.5, ymin = -0.90000000000000002, ymax = -0.80000000000000004), fill = "#75C2F6", alpha = 0.29999999999999999) + 
        ggplot2::geom_text(data = beta_sheets, ggplot2::aes(x = (xstart + xend)/2, y = -0.84999999999999998, label = col), 
            size = 8/points_per_mm, color = "#A31300") + ggplot2::geom_rect(data = alpha_helices, ggplot2::aes(xmin = xstart - 
        0.5, xmax = xend + 0.5, ymin = -0.90000000000000002, ymax = -0.80000000000000004), fill = "#C68EFD", alpha = 0.29999999999999999) + 
        ggplot2::geom_text(data = alpha_helices, ggplot2::aes(x = (xstart + xend)/2, y = -0.84999999999999998, label = col), 
            size = 8/points_per_mm, color = "#A31300") + ggplot2::scale_x_continuous(breaks = seq(0, sequence_length, by = 10), 
        expand = c(0.02, 0.02)) + ggplot2::scale_y_continuous(limits = c(-1.2, 0.5)) + ggplot2::labs(title = "KRAS Protein Sequence Annotation", 
        x = "Amino Acid Position", y = "", caption = "Visualization of KRAS structural and functional features") + ggplot2::theme_minimal(base_size = 8) + 
        ggplot2::theme(text = ggplot2::element_text(size = 8), axis.text.x = ggplot2::element_text(size = 8), axis.text.y = ggplot2::element_blank(), 
            axis.ticks.y = ggplot2::element_blank(), axis.title = ggplot2::element_text(size = 8), plot.title = ggplot2::element_text(size = 8), 
            plot.caption = ggplot2::element_text(size = 8, hjust = 0.5), legend.title = ggplot2::element_text(size = 8), 
            legend.text = ggplot2::element_text(size = 8), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
            panel.background = ggplot2::element_blank(), legend.position = "bottom")
    return(p)
}

