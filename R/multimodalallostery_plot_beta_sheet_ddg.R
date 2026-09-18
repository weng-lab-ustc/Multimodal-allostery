#' plot beta sheet ddg
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_beta_sheet_ddg()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param data Input or option used by this workflow; see Usage for the exact interface.
#' @param y_label Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_beta_sheet_ddg <- function (data, y_label) 
{
    ggplot2::ggplot(data,ggplot2::aes(x = colors_type, y = ddG)) + ggplot2::geom_violin() + ggplot2::geom_jitter(size = 0.34999999999999998, 
        height = 0) + ggplot2::ylab(y_label) + ggplot2::xlab("beta sheet") + ggplot2::theme_classic() + ggplot2::theme(text = ggplot2::element_text(size = 10), 
        axis.text = ggplot2::element_text(size = 10), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.text = ggplot2::element_text(size = 10))
}

