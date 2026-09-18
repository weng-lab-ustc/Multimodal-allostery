#' .ma aes string
#'
#' Workflow-only September 2026 implementation of `.ma_aes_string()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param ... Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
.ma_aes_string <- function (...) 
{
    expressions <- list(...)
    mapping <- ggplot2::aes()
    for (name in names(expressions)) {
        value <- expressions[[name]]
        canonical_name <- if (name == "color") 
            "colour"
        else name
        if (!is.null(value)) 
            mapping[[canonical_name]] <- rlang::new_quosure(rlang::parse_expr(value), env = parent.frame())
    }
    mapping
}

