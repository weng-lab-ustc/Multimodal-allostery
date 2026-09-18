#' add region marker
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_add_region_marker()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param data Input or option used by this workflow; see Usage for the exact interface.
#' @param nbp_residues Input or option used by this workflow; see Usage for the exact interface.
#' @param switch_ii_residues Input or option used by this workflow; see Usage for the exact interface.
#' @param switch_i_residues Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_add_region_marker <- function (data, nbp_residues = c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 
    146, 147), switch_ii_residues = c(58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76), switch_i_residues = c(25, 
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40)) 
{
    data[, `:=`(region, "Other")]
    data[Pos_real %in% nbp_residues, `:=`(region, "NBP")]
    data[region != "NBP" & Pos_real %in% switch_i_residues, `:=`(region, "Switch I")]
    data[region != "NBP" & Pos_real %in% switch_ii_residues, `:=`(region, "Switch II")]
    anticorrelated_types <- c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y")
    data[, `:=`(plot_group, as.character(final_classification))]
    data[final_classification %in% anticorrelated_types & region == "NBP", `:=`(plot_group, "Anticorrelated_NBP")]
    data[final_classification %in% anticorrelated_types & region == "Switch I", `:=`(plot_group, "Anticorrelated_SwitchI")]
    data[final_classification %in% anticorrelated_types & region == "Switch II", `:=`(plot_group, "Anticorrelated_SwitchII")]
    data[final_classification %in% anticorrelated_types & region == "Other", `:=`(plot_group, "Anticorrelated_Other")]
    return(data)
}

