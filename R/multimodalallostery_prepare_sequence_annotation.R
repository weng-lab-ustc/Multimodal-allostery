#' prepare sequence annotation
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_prepare_sequence_annotation()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param sequence Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_sites_1 Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_sites_2 Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_color_1 Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_color_2 Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_prepare_sequence_annotation <- function (sequence, binding_sites_1, binding_sites_2, binding_color_1 = "#1B38A6", binding_color_2 = "#F4270C") 
{
    sequence_length <- nchar(sequence)
    positions <- seq_len(sequence_length)
    result <- data.frame(position = positions, residue = strsplit(sequence, "", fixed = TRUE)[[1]], binding_1 = positions %in% 
        binding_sites_1, binding_2 = positions %in% binding_sites_2)
    result$binding_group <- interaction(result$binding_1, result$binding_2)
    return(list(sequence_data = result, sequence_length = sequence_length, binding_colors = c("black", binding_color_1, binding_color_2, 
        "#F4AD0C")))
}

