#' identify mutation positions
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_identify_mutation_positions()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input Input or option used by this workflow; see Usage for the exact interface.
#' @param wt_aa Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_identify_mutation_positions <- function (input, wt_aa) 
{
    output <- input
    output[, `:=`(AA_Pos1, which(unlist(strsplit(aa_seq, "")) != unlist(strsplit(wt_aa, ""))[1:nchar(aa_seq)])[1]), aa_seq]
    output[, `:=`(AA_Pos2, which(unlist(strsplit(aa_seq, "")) != unlist(strsplit(wt_aa, ""))[1:nchar(aa_seq)])[2]), aa_seq]
    for (i in 1:188) {
        output[AA_Pos1 == i, `:=`(mt1, substr(aa_seq, i, i))]
    }
    for (i in 1:188) {
        output[AA_Pos2 == i, `:=`(mt2, substr(aa_seq, i, i))]
    }
    for (i in 1:188) {
        output[AA_Pos1 == i, `:=`(wtcodon1, substr(wt_aa, i, i))]
    }
    for (i in 1:188) {
        output[AA_Pos2 == i, `:=`(wtcodon2, substr(wt_aa, i, i))]
    }
    output[, `:=`(codon1, substr(aa_seq, AA_Pos1, AA_Pos1))]
    output[, `:=`(codon2, substr(aa_seq, AA_Pos2, AA_Pos2))]
    output[, `:=`(AA_Pos1, AA_Pos1 + 1)]
    output[, `:=`(AA_Pos2, AA_Pos2 + 1)]
    output[, `:=`(mt1, paste0(wtcodon1, AA_Pos1, codon1))]
    output[, `:=`(mt2, paste0(wtcodon2, AA_Pos2, codon2))]
    return(output)
}

