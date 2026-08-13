#' Load Dd g Data.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input Input file path or input data object.
#' @param assay_sele Selected assay identifier.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_load_dd_g_data <- function(input, assay_sele) {
    ddG <- fread(input)
    ddG[, `:=`(Pos_real, Pos_ref + 1)]
    ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
    ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
    ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
    ddG <- ddG[id != "WT"]
    result <- ddG[, .(mt, Pos_real, `mean_kcal/mol`)]
    setnames(result, "mean_kcal/mol", "ddG")
    result[, `:=`(assay, assay_sele)]
    return(result)
}

