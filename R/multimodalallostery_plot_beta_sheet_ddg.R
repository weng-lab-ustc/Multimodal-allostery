#' Plot Beta Sheet Ddg.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data Input data frame or data table.
#' @param y_label Value supplied for `y_label`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_beta_sheet_ddg <- function(data, y_label) {
    ggplot2::ggplot(data, ggplot2::aes(x = colors_type, y = `mean_kcal/mol`)) + ggplot2::geom_violin() + ggplot2::geom_jitter(size = 0.35, 
        height = 0) + ggplot2::ylab(y_label) + ggplot2::xlab("beta sheet") + ggplot2::theme_classic() + ggplot2::theme(text = ggplot2::element_text(size = 10), 
        axis.text = ggplot2::element_text(size = 10), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.text = ggplot2::element_text(size = 10))
}

