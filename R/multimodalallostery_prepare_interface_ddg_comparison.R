#' prepare interface ddg comparison
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_prepare_interface_ddg_comparison()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param data_x Input or option used by this workflow; see Usage for the exact interface.
#' @param data_y Input or option used by this workflow; see Usage for the exact interface.
#' @param interface_x Input or option used by this workflow; see Usage for the exact interface.
#' @param interface_y Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_x Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_y Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_prepare_interface_ddg_comparison <- function(
    data_x, data_y, interface_x, interface_y, assay_x, assay_y
) {x <- data_x[Pos_real %in% interface_x, .(mt, Pos_real, ddG, `std_kcal/mol`)]
   y <- data_y[Pos_real %in% interface_y, .(mt, Pos_real, ddG, `std_kcal/mol`)]
   data.table::setnames(x, c("ddG", "std_kcal/mol"), c(assay_x, "std_kcal/mol.x"))
   data.table::setnames(y, c("ddG", "std_kcal/mol"), c(assay_y, "std_kcal/mol.y"))
   result <- merge(x, y, by = c("mt", "Pos_real"), all = FALSE)
   if (!all(c(assay_x, assay_y) %in% names(result)))
        stop("Prepared data do not contain the requested assay columns.")
   result
}

