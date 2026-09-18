#' Plot Ddg Fitness per Block.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param pre_nor Value supplied for `pre_nor`.
#' @param phenotypen Value supplied for `phenotypen`.
#' @param rotate_x_axis Value supplied for `rotate_x_axis`.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_plot_ddg_fitness_per_block <- function(pre_nor = pre_nor, phenotypen = phenotypen, rotate_x_axis = TRUE) {
    pre_nor
    lm_mochi <- lm(pre_nor_fitness ~ ob_nor_fitness, pre_nor[phenotype == phenotypen, ])
    p <- ggplot2::ggplot() + ggplot2::stat_binhex(data = pre_nor[phenotype == phenotypen, ], ggplot2::aes(x = ob_nor_fitness, 
        y = pre_nor_fitness), bins = 50, size = 0, color = "black") + ggplot2::scale_fill_gradient(low = "white", high = "black", 
        trans = "log10", guide = ggplot2::guide_colorbar(barwidth = 0.5, barheight = 1.5)) + ggplot2::geom_hline(yintercept = 0) + 
        ggplot2::geom_vline(xintercept = 0) + ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") + ggplot2::annotate("text", 
        x = -0.8, y = 0.3, label = paste0("R\u00B2 = ", round(summary(lm_mochi)$r.squared, 2)), size = 8 * 0.35) + ggplot2::theme_classic() + 
        ggplot2::xlab("Observed fitness") + ggplot2::ylab("Predicted fitness") + ggplot2::theme(text = ggplot2::element_text(size = 8), 
        axis.text = ggplot2::element_text(size = 8), legend.text = ggplot2::element_text(size = 8), legend.key.size = ggplot2::unit(1, 
            "cm"), plot.title = ggplot2::element_text(size = 8)) + ggplot2::coord_fixed()
    if (rotate_x_axis) {
        p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
    }
    return(p)
}

