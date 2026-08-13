#' Classify Mutation Type.
#'
#' Reusable classification/annotation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data Input data frame or data table.
#'
#' @return The input data with classifications or annotations added.
#' @export
multimodalallostery_classify_mutation_type <- function(data) {
    classified_data <- data %>% dplyr::mutate(mut_type = dplyr::case_when(Nham_aa == 0 & Nham_nt > 0 ~ "Synonymous", STOP == 
        TRUE | STOP_readthrough == TRUE ~ "Stop", Nham_aa > 0 & indel == FALSE & STOP == FALSE & STOP_readthrough == FALSE ~ 
        "Missense")) %>% dplyr::filter(!is.na(mut_type)) %>% dplyr::mutate(mut_type = factor(mut_type, levels = c("Synonymous", 
        "Missense", "Stop")))
    return(classified_data)
}

