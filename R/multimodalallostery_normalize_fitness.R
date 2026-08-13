#' Normalize Fitness.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param block1 Value supplied for `block1`.
#' @param block2 Value supplied for `block2`.
#' @param block3 Value supplied for `block3`.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_normalize_fitness <- function(block1, block2, block3) {
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

