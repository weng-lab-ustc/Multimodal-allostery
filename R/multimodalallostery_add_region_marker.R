#' Add Region Marker.
#'
#' Reusable utility logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data Input data frame or data table.
#' @param nbp_residues Residue positions assigned to the nucleotide-binding pocket.
#' @param switch_i_residues Residue positions assigned to Switch I.
#' @param switch_ii_residues Residue positions assigned to Switch II.
#'
#' @return The transformed result produced by the function.
#' @export
multimodalallostery_add_region_marker <- function(data, nbp_residues, switch_i_residues, switch_ii_residues) {
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

