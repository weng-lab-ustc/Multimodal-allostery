#' Identify K13 versus K19 mutations that promote both assays
#'
#' @param input_k13 Path to the K13 input file.
#' @param input_k19 Path to the K19 input file.
#' @param anno Annotation data or annotation file path.
#' @param binding_sites_map Named list mapping assays to binding-site residues.
#' @param nbp_residues Residue positions in the nucleotide-binding pocket.
#' @return A data table containing mutations classified as promoting both assays.
#' @export
multimodalallostery_get_k13_vs_k19_both_promoting <- function(input_k13, input_k19, anno,
                                          binding_sites_map, nbp_residues) {
  prepared <- multimodalallostery_prepare_mapped_merged_data_with_fdr(
    input_k13, input_k19, "K13", "K19", anno, binding_sites_map
  )
  merged_data <- prepared$data
  merged_data[, pass_FDR_K13 := p_adj_x < 0.05]
  merged_data[, pass_FDR_K19 := p_adj_y < 0.05]
  merged_data[, direction_class := multimodalallostery_classify_by_direction(
    ddG_K13, ddG_K19, prepared$threshold_x, prepared$threshold_y
  )]
  merged_data[, final_classification := multimodalallostery_reclassify_by_fdr(
    direction_class, pass_FDR_K13, pass_FDR_K19
  )]
  result <- merged_data[final_classification == "Both promoting"]
  result[, is_NBP := Pos_real %in% nbp_residues]
  result
}
