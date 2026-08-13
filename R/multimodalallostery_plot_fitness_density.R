#' Plot Fitness Density.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data Input data frame or data table.
#' @param assay_type Value supplied for `assay_type`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_fitness_density <- function(data, assay_type) {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = fitness_normalized, color = mut_type)) + ggplot2::geom_density(linewidth = 1) + 
        ggplot2::scale_color_manual(values = c(Synonymous = "#09B636", Missense = "#F4AD0C", Stop = "#FF6A56")) + ggplot2::labs(title = paste0(toupper(assay_type), 
        " - Fitness Distribution"), x = "Normalized Fitness", y = "Density", color = "Mutation Type") + ggplot2::xlim(-1.5, 
        0.5) + ggplot2::theme_classic() + ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5, 
        size = 12), text = ggplot2::element_text(size = 10), legend.title = ggplot2::element_text(size = 10), legend.text = ggplot2::element_text(size = 9))
    return(p)
}

