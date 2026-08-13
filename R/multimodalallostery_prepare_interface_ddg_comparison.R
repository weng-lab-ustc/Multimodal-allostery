#' Prepare Interface Ddg Comparison.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data_x Value supplied for `data_x`.
#' @param data_y Value supplied for `data_y`.
#' @param interface_x Value supplied for `interface_x`.
#' @param interface_y Value supplied for `interface_y`.
#' @param assay_x Identifier for the first assay.
#' @param assay_y Identifier for the second assay.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_interface_ddg_comparison <- function(data_x, data_y, interface_x, interface_y, assay_x, assay_y) {
    x <- data_x[Pos_real %in% interface_x, c(1:3, 23:27)] %>% tidyr::spread(key = assay, value = `mean_kcal/mol`)
    y <- data_y[Pos_real %in% interface_y, c(1:3, 23:27)] %>% tidyr::spread(key = assay, value = `mean_kcal/mol`)
    result <- merge(x, y, by = c("mt", "Pos_real"), all = FALSE)
    required <- c(assay_x, assay_y)
    if (!all(required %in% names(result))) {
        stop("Prepared data do not contain the requested assay columns.")
    }
    return(result)
}

