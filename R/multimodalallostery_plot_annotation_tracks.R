#' Plot Annotation Tracks.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data Input data frame or data table.
#' @param residue_levels Value supplied for `residue_levels`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_annotation_tracks <- function(data, residue_levels) {
    ggplot2::ggplot(data, ggplot2::aes(x = residue_order, y = track)) + ggplot2::geom_tile(width = 0.95, height = 0.75, fill = "grey70") + 
        ggplot2::scale_x_discrete(drop = FALSE, breaks = residue_levels, labels = residue_levels) + ggplot2::labs(x = "KRAS residue position (ordered by Pearson's R)") + 
        ggplot2::theme_classic(base_size = 11) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, 
        vjust = 0.5, size = 6), axis.ticks.x = ggplot2::element_line(linewidth = 0.25), axis.text.y = ggplot2::element_text(size = 9), 
        plot.margin = ggplot2::margin(0, 5, 5, 5))
}

