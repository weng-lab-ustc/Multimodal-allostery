#' load site mutation data
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_load_site_mutation_data()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_sele Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_load_site_mutation_data <- function (input, assay_sele) 
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
    mut_data <- ddG[id != "WT", .(mt = paste0(wt_codon, Pos_real, mt_codon), Pos_real, wt_codon, mt_codon, ddG = `mean_kcal/mol`, 
        ddG_std = `std_kcal/mol`)]
    all_positions <- unique(mut_data$Pos_real)
    wt_rows <- data.table::data.table(mt = paste0("WT", all_positions, "WT"), Pos_real = all_positions, wt_codon = "WT", 
        mt_codon = "WT", ddG = 0, ddG_std = 0)
    result <- rbind(mut_data, wt_rows)
    result[, `:=`(assay, assay_sele)]
    data.table::setorder(result, Pos_real, mt_codon)
    return(result)
}

