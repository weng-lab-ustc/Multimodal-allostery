#' load mapped mutation data
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_load_mapped_mutation_data()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_sele Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_load_mapped_mutation_data <- function (input, assay_sele) 
{
    if (is.character(input)) {
        ddG <- .ma_read_table(input)
    }
    else {
        ddG <- data.table::as.data.table(input)
    }
    ddG[, `:=`(Pos_real, Pos_ref + 1)]
    ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
    ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
    ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
    ddG <- ddG[id != "WT"]
    result <- ddG[, .(mt, Pos_real, wt_codon, mt_codon, `mean_kcal/mol`, `std_kcal/mol`)]
    data.table::setnames(result, "mean_kcal/mol", "ddG")
    data.table::setnames(result, "std_kcal/mol", "ddG_std")
    result[, `:=`(assay, assay_sele)]
    result[, `:=`(mt, toupper(trimws(mt)))]
    return(result)
}

