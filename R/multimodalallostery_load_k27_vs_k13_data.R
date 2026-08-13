#' Load and prepare K27 versus K13 data
#'
#' @param input_k27 Path to the K27 input file.
#' @param input_k13 Path to the K13 input file.
#' @param anno Annotation data or annotation file path.
#' @param binding_sites_map Named list mapping assays to binding-site residues.
#' @return A list containing merged data and the K27/K13 thresholds.
#' @export
multimodalallostery_load_k27_vs_k13_data <- function(input_k27, input_k13, anno, binding_sites_map) {
  prepared <- multimodalallostery_prepare_mapped_merged_data_with_fdr(
    input_k27, input_k13, "K27", "K13", anno, binding_sites_map
  )
  list(data = prepared$data, threshold_K27 = prepared$threshold_x, threshold_K13 = prepared$threshold_y)
}
