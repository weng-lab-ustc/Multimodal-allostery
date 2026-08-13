#' Plot Sequence Annotation.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param prepared_sequence Value supplied for `prepared_sequence`.
#' @param gtp_pocket Value supplied for `gtp_pocket`.
#' @param functional_loops Value supplied for `functional_loops`.
#' @param core_residues Value supplied for `core_residues`.
#' @param beta_sheets Value supplied for `beta_sheets`.
#' @param alpha_helices Value supplied for `alpha_helices`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_sequence_annotation <- function(prepared_sequence, gtp_pocket, functional_loops, core_residues, beta_sheets, alpha_helices) {
    sequence_data <- prepared_sequence$sequence_data
    sequence_length <- prepared_sequence$sequence_length
    points_per_mm <- 72.27/25.4
    p <- ggplot2::ggplot() + ggplot2::geom_text(data = sequence_data, ggplot2::aes(x = position, y = 0, label = residue, 
        color = binding_group), size = 8/points_per_mm, family = "mono", angle = 90, vjust = 0.5, hjust = 0.5) + ggplot2::scale_color_manual(values = prepared_sequence$binding_colors, 
        labels = c("Other", "K13 binding", "RAF1 binding", "Both"), name = "Binding Sites") + ggplot2::geom_rect(data = data.frame(position = gtp_pocket), 
        ggplot2::aes(xmin = position - 0.5, xmax = position + 0.5, ymin = -0.3, ymax = -0.2), fill = "#F1DD10", alpha = 0.5) + 
        ggplot2::annotate("text", x = mean(range(gtp_pocket)), y = -0.25, label = "GTP Pocket", color = "#F1DD10", size = 8/points_per_mm) + 
        ggplot2::geom_rect(data = functional_loops, ggplot2::aes(xmin = xstart - 0.5, xmax = xend + 0.5, ymin = -0.5, ymax = -0.4, 
            fill = col), alpha = 0.3) + ggplot2::geom_text(data = functional_loops, ggplot2::aes(x = (xstart + xend)/2, y = -0.45, 
        label = col), size = 8/points_per_mm, color = "black") + ggplot2::scale_fill_manual(values = c("#FFB0A5", "#FF0066", 
        "#007A20"), name = "Functional Loops") + ggplot2::geom_rect(data = data.frame(position = core_residues), ggplot2::aes(xmin = position - 
        0.5, xmax = position + 0.5, ymin = -0.7, ymax = -0.6), fill = "#85409D", alpha = 0.9) + ggplot2::annotate("text", 
        x = mean(range(core_residues)), y = -0.65, label = "Core", color = "black", size = 8/points_per_mm) + ggplot2::geom_rect(data = beta_sheets, 
        ggplot2::aes(xmin = xstart - 0.5, xmax = xend + 0.5, ymin = -0.9, ymax = -0.8), fill = "#75C2F6", alpha = 0.3) + 
        ggplot2::geom_text(data = beta_sheets, ggplot2::aes(x = (xstart + xend)/2, y = -0.85, label = col), size = 8/points_per_mm, 
            color = "#A31300") + ggplot2::geom_rect(data = alpha_helices, ggplot2::aes(xmin = xstart - 0.5, xmax = xend + 
        0.5, ymin = -0.9, ymax = -0.8), fill = "#C68EFD", alpha = 0.3) + ggplot2::geom_text(data = alpha_helices, ggplot2::aes(x = (xstart + 
        xend)/2, y = -0.85, label = col), size = 8/points_per_mm, color = "#A31300") + ggplot2::scale_x_continuous(breaks = seq(0, 
        sequence_length, by = 10), expand = c(0.02, 0.02)) + ggplot2::scale_y_continuous(limits = c(-1.2, 0.5)) + ggplot2::labs(title = "KRAS Protein Sequence Annotation", 
        x = "Amino Acid Position", y = "", caption = "Visualization of KRAS structural and functional features") + ggplot2::theme_minimal(base_size = 8) + 
        ggplot2::theme(text = ggplot2::element_text(size = 8), axis.text.x = ggplot2::element_text(size = 8), axis.text.y = ggplot2::element_blank(), 
            axis.ticks.y = ggplot2::element_blank(), axis.title = ggplot2::element_text(size = 8), plot.title = ggplot2::element_text(size = 8), 
            plot.caption = ggplot2::element_text(size = 8, hjust = 0.5), legend.title = ggplot2::element_text(size = 8), 
            legend.text = ggplot2::element_text(size = 8), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
            panel.background = ggplot2::element_blank(), legend.position = "bottom")
    return(p)
}

