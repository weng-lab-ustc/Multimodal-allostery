#' Plot Residue Correlations.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data Input data frame or data table.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_residue_correlations <- function(data) {
    ggplot2::ggplot(data, ggplot2::aes(x = residue_order, y = R)) + ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, 
        colour = "grey40") + ggplot2::geom_point(ggplot2::aes(fill = sig_group), shape = 21, size = 2.3, colour = "white", 
        stroke = 0.25) + ggplot2::scale_fill_manual(values = c(Positive = "#F4AD0C", Negative = "#75C2F6", NS = "grey80")) + 
        ggplot2::scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5), expand = ggplot2::expansion(mult = c(0.02, 
            0.02))) + ggplot2::labs(y = "Residue Correlation Pearson's R between RAF1 and K13", fill = NULL) + ggplot2::theme_classic(base_size = 11) + 
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "top", 
            plot.margin = ggplot2::margin(5, 5, 0, 5))
}

