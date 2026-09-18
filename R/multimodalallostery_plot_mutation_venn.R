#' plot mutation venn
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_mutation_venn()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param mutation_sets Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_mutation_venn <- function (mutation_sets) 
{
    ggVennDiagram::ggVennDiagram(mutation_sets, label_alpha = 0, edge_size = 0, label_size = 5, set_size = 5) + ggplot2::scale_fill_gradient(low = "white", 
        high = "#75C2F6") + ggplot2::scale_color_manual(values = rep("transparent", length(mutation_sets))) + ggplot2::theme_void() + 
        ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 12))
}

