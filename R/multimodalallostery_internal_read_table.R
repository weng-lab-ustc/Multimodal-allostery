#' .ma read table
#'
#' Workflow-only September 2026 implementation of `.ma_read_table()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input Input or option used by this workflow; see Usage for the exact interface.
#' @param ... Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
.ma_read_table <- function (input, ...) 
{
    if (is.character(input)) 
        data.table::fread(input, ...)
    else data.table::copy(data.table::as.data.table(input))
}

