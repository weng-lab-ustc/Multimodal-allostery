#' prepare beta sheet ddg
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_prepare_beta_sheet_ddg()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param ddg Input or option used by this workflow; see Usage for the exact interface.
#' @param annotation Input or option used by this workflow; see Usage for the exact interface.
#' @param beta_sheet_ranges Input or option used by this workflow; see Usage for the exact interface.
#' @param sheet_levels Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_prepare_beta_sheet_ddg <- function (ddg, annotation, beta_sheet_ranges, sheet_levels = c("b2", "b3", "b1", "b4", "b5", "b6")) 
{
    ddg <- data.table::copy(data.table::as.data.table(ddg))
    ddg[, Pos := Pos_real]
    result <- merge(ddg, annotation, by = "Pos", all = TRUE)
    result[, `:=`(colors_type, "others")]
    for (sheet in beta_sheet_ranges$col) {
        bounds <- beta_sheet_ranges[beta_sheet_ranges$col == sheet, ]
        result[Pos >= bounds$xstart & Pos <= bounds$xend, `:=`(colors_type, sheet)]
    }
    result <- result[colors_type != "others"]
    result[, `:=`(colors_type, factor(colors_type, levels = sheet_levels))]
    result
}

