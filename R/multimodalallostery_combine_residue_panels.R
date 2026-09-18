#' combine residue panels
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_combine_residue_panels()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param correlation_plot Input or option used by this workflow; see Usage for the exact interface.
#' @param annotation_plot Input or option used by this workflow; see Usage for the exact interface.
#' @param heights Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_combine_residue_panels <- function (correlation_plot, annotation_plot, heights = c(4, 2.5)) 
{
    patchwork::wrap_plots(correlation_plot, annotation_plot, ncol = 1, heights = heights)
}

