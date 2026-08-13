#' Prepare Beta Sheet Ddg.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param ddg Value supplied for `ddg`.
#' @param annotation Annotation data frame.
#' @param beta_sheet_ranges Value supplied for `beta_sheet_ranges`.
#' @param sheet_levels Value supplied for `sheet_levels`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_beta_sheet_ddg <- function(ddg, annotation, beta_sheet_ranges, sheet_levels = c("b2", "b3", "b1", "b4", "b5", "b6")) {
    ddg <- data.table::copy(data.table::as.data.table(ddg))
    ddg[, `:=`(Pos, Pos_ref + 1)]
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

