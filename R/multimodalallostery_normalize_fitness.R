#' normalize fitness
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_normalize_fitness()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param block1 Input or option used by this workflow; see Usage for the exact interface.
#' @param block2 Input or option used by this workflow; see Usage for the exact interface.
#' @param block3 Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_normalize_fitness <- function (block1, block2, block3) 
{
    nor_fit <- wlab.block::nor_fitness(block1 = block1, block2 = block2, block3 = block3)
    WT_mean <- nor_fit %>% dplyr::filter(WT == TRUE) %>% dplyr::summarise(WT = sum(nor_fitness/nor_fitness_sigma^2, na.rm = TRUE)/sum(1/nor_fitness_sigma^2, 
        na.rm = TRUE)) %>% dplyr::pull(WT)
    STOP_mean <- nor_fit %>% dplyr::filter(STOP == TRUE) %>% dplyr::summarise(STOP = sum(nor_fitness/nor_fitness_sigma^2, 
        na.rm = TRUE)/sum(1/nor_fitness_sigma^2, na.rm = TRUE)) %>% dplyr::pull(STOP)
    cat("WT reference =", WT_mean, "\n")
    cat("STOP reference =", STOP_mean, "\n")
    nor_fit <- nor_fit %>% dplyr::mutate(fitness_normalized = (nor_fitness - WT_mean)/(WT_mean - STOP_mean))
    return(nor_fit)
}

