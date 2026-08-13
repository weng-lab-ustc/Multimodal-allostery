#' Plot Allosteric Pair Counts.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data Input data frame or data table.
#' @param colors Value supplied for `colors`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_allosteric_pair_counts <- function(data, colors) {
    ggplot2::ggplot(data, ggplot2::aes(x = Pair, y = Count, fill = Type)) + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(0.7), 
        width = 0.6, color = "white", linewidth = 0.3) + ggplot2::scale_fill_manual(values = colors) + ggplot2::labs(title = "Allosteric Mutations Exclusive to Each Binder Pair", 
        x = "Binder Pair", y = "Number of Allosteric Mutations", fill = "Type") + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
        hjust = 1, size = 10), axis.text.y = ggplot2::element_text(size = 10), axis.title = ggplot2::element_text(size = 12), 
        axis.line = ggplot2::element_line(color = "black", linewidth = 0.5), axis.ticks = ggplot2::element_line(color = "black", 
            linewidth = 0.5), plot.title = ggplot2::element_text(size = 14, hjust = 0.5), legend.position = "top", legend.title = ggplot2::element_text(size = 11), 
        legend.text = ggplot2::element_text(size = 10), panel.background = ggplot2::element_rect(fill = "white"), plot.background = ggplot2::element_rect(fill = "white"))
}

