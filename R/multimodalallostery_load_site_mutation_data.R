#' Load Site Mutation Data.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input Input file path or input data object.
#' @param assay_sele Selected assay identifier.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_load_site_mutation_data <- function(input, assay_sele) {
    if (is.character(input)) {
        ddG <- fread(input)
    }
    else {
        ddG <- as.data.table(input)
    }
    ddG[, `:=`(Pos_real, Pos_ref + 1)]
    ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
    ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
    mut_data <- ddG[id != "WT", .(mt = paste0(wt_codon, Pos_real, mt_codon), Pos_real, wt_codon, mt_codon, ddG = `mean_kcal/mol`, 
        ddG_std = `std_kcal/mol`)]
    all_positions <- unique(mut_data$Pos_real)
    wt_rows <- data.table(mt = paste0("WT", all_positions, "WT"), Pos_real = all_positions, wt_codon = "WT", mt_codon = "WT", 
        ddG = 0, ddG_std = 0)
    result <- rbind(mut_data, wt_rows)
    result[, `:=`(assay, assay_sele)]
    setorder(result, Pos_real, mt_codon)
    return(result)
}

