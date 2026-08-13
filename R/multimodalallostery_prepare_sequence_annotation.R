#' Prepare Sequence Annotation.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param sequence Value supplied for `sequence`.
#' @param binding_sites_1 Value supplied for `binding_sites_1`.
#' @param binding_sites_2 Value supplied for `binding_sites_2`.
#' @param binding_color_1 Value supplied for `binding_color_1`.
#' @param binding_color_2 Value supplied for `binding_color_2`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_sequence_annotation <- function(sequence, binding_sites_1, binding_sites_2, binding_color_1 = "#1B38A6", binding_color_2 = "#F4270C") {
    sequence_length <- nchar(sequence)
    positions <- seq_len(sequence_length)
    result <- data.frame(position = positions, residue = strsplit(sequence, "", fixed = TRUE)[[1]], binding_1 = positions %in% 
        binding_sites_1, binding_2 = positions %in% binding_sites_2)
    result$binding_group <- interaction(result$binding_1, result$binding_2)
    return(list(sequence_data = result, sequence_length = sequence_length, binding_colors = c("black", binding_color_1, binding_color_2, 
        "#F4AD0C")))
}

