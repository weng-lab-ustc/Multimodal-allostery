#' .ma save plot
#'
#' Workflow-only September 2026 implementation of `.ma_save_plot()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param filename Input or option used by this workflow; see Usage for the exact interface.
#' @param ... Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
.ma_save_plot <- function (filename, ...) 
{
    if (file.exists(filename)) 
        stop("Refusing to overwrite existing output: ", filename)
    ggplot2::ggsave(filename = filename, ...)
}

