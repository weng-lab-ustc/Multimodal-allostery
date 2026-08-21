#' Plot Fitness Correlation Blocks.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input Input file path or input data object.
#' @param assay_name Value supplied for `assay_name`.
#' @param colour_scheme Value supplied for `colour_scheme`.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_plot_fitness_correlation_blocks_BI1 <- function(input, assay_name, colour_scheme) {
    d <- GGally::ggpairs(input, columns = c(18, 20, 22), columnLabels = c("replicate 1", "replicate 2", "replicate 3"), mapping = ggplot2::aes(color = block), 
        lower = list(continuous = function(data, mapping, ...) {
            ggplot2::ggplot(data = data, mapping = mapping) + ggplot2::geom_bin2d(bins = 100, alpha = 0.2) + ggplot2::scale_fill_gradient(low = "white", 
                high = "black")
        }), upper = list(continuous = GGally::wrap("cor", mapping = ggplot2::aes(color = block), size = 8 * 0.35)), diag = list(continuous = "blankDiag")) + 
        ggplot2::scale_color_manual(values = c(colour_scheme[["red"]], colour_scheme[["green"]], colour_scheme[["blue"]])) + 
        ggplot2::ggtitle(assay_name) + ggpubr::theme_classic2() + ggplot2::theme(text = ggplot2::element_text(size = 8), 
        axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = ggplot2::element_text(size = 8), 
        axis.title = ggplot2::element_text(size = 8), strip.text.x = ggplot2::element_text(size = 8), strip.text.y = ggplot2::element_text(size = 8), 
        legend.text = ggplot2::element_text(size = 8), plot.title = ggplot2::element_text(hjust = 0.5, size = 8))
    return(d)
}

