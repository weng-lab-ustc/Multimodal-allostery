#' plot allosteric pair counts
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_allosteric_pair_counts()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param data Input or option used by this workflow; see Usage for the exact interface.
#' @param colors Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_allosteric_pair_counts <- function (data, colors) 
{
    ggplot2::ggplot(data, ggplot2::aes(x = Pair, y = Count, fill = Type)) + ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(0.69999999999999996), 
        width = 0.59999999999999998, color = "white", linewidth = 0.29999999999999999) + ggplot2::scale_fill_manual(values = colors) + 
        ggplot2::labs(title = "Allosteric Mutations Exclusive to Each Binder Pair", x = "Binder Pair", y = "Number of Allosteric Mutations", 
            fill = "Type") + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, 
        size = 10), axis.text.y = ggplot2::element_text(size = 10), axis.title = ggplot2::element_text(size = 12), axis.line = ggplot2::element_line(color = "black", 
        linewidth = 0.5), axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.5), plot.title = ggplot2::element_text(size = 14, 
        hjust = 0.5), legend.position = "top", legend.title = ggplot2::element_text(size = 11), legend.text = ggplot2::element_text(size = 10), 
        panel.background = ggplot2::element_rect(fill = "white"), plot.background = ggplot2::element_rect(fill = "white"))
}

