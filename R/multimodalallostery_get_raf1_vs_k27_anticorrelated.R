#' Identify RAF1 versus K27 anticorrelated mutations
#'
#' @param input_raf1 Path to the RAF1 input file.
#' @param input_k27 Path to the K27 input file.
#' @param anno Annotation data or annotation file path.
#' @param binding_sites_map Named list mapping assays to binding-site residues.
#' @param nbp_residues Residue positions in the nucleotide-binding pocket.
#' @return A data table containing anticorrelated mutations.
#' @export
multimodalallostery_get_raf1_vs_k27_anticorrelated <- function(input_raf1, input_k27, anno,
                                           binding_sites_map, nbp_residues) {
  prepared <- multimodalallostery_prepare_mapped_merged_data_with_fdr(
    input_raf1, input_k27, "RAF1", "K27", anno, binding_sites_map
  )
  merged_data <- prepared$data
  merged_data[, pass_FDR_RAF1 := p_adj_x < 0.05]
  merged_data[, pass_FDR_K27 := p_adj_y < 0.05]
  merged_data[, direction_class := multimodalallostery_classify_by_direction(
    ddG_RAF1, ddG_K27, prepared$threshold_x, prepared$threshold_y
  )]
  merged_data[, final_classification := multimodalallostery_reclassify_by_fdr(
    direction_class, pass_FDR_RAF1, pass_FDR_K27
  )]
  result <- merged_data[final_classification %in% c(
    "Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y"
  )]
  result[, is_NBP := Pos_real %in% nbp_residues]
  result
}
