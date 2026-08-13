#' Plot Mutation Venn.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param mutation_sets Value supplied for `mutation_sets`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_mutation_venn <- function(mutation_sets) {
    ggVennDiagram::ggVennDiagram(mutation_sets, label_alpha = 0, edge_size = 0, label_size = 5, set_size = 5) + ggplot2::scale_fill_gradient(low = "white", 
        high = "#75C2F6") + ggplot2::scale_color_manual(values = rep("transparent", length(mutation_sets))) + ggplot2::theme_void() + 
        ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 12))
}

