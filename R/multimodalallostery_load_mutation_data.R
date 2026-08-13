#' Load Mutation Data.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input Input file path or input data object.
#' @param assay_sele Selected assay identifier.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_load_mutation_data <- function(input, assay_sele) {
    if (is.character(input)) {
        ddG <- fread(input)
    }
    else {
        ddG <- as.data.table(input)
    }
    ddG[, `:=`(Pos_real, Pos_ref + 1)]
    ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
    ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
    ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
    ddG <- ddG[id != "WT"]
    result <- ddG[, .(mt, Pos_real, wt_codon, mt_codon, `mean_kcal/mol`, `std_kcal/mol`)]
    setnames(result, "mean_kcal/mol", "ddG")
    setnames(result, "std_kcal/mol", "ddG_std")
    result[, `:=`(assay, assay_sele)]
    return(result)
}

