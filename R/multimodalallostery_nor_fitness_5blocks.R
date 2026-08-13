#' Nor Fitness 5blocks.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param required_file Value supplied for `required_file`.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_nor_fitness_5blocks <- function(required_file) {
    if (length(required_file) != 5) 
        stop("Five file paths need to be provided.")
    all_data <- list()
    gr_lm <- data.frame(matrix(nrow = 2, ncol = 0))
    fit_lm <- data.frame(matrix(nrow = 2, ncol = 0))
    for (i in 1:5) {
        load(required_file[i])
        all_variants$gr_over_sigmasquared <- all_variants$growthrate/(all_variants$growthrate_sigma)^2
        all_variants$one_over_sigmasquared <- 1/(all_variants$growthrate_sigma)^2
        stop1 <- sum(all_variants[STOP == TRUE, ]$gr_over_sigmasquared, na.rm = TRUE)/sum(all_variants[STOP == TRUE, ]$one_over_sigmasquared, 
            na.rm = TRUE)
        wt1 <- sum(all_variants[WT == TRUE, ]$gr_over_sigmasquared, na.rm = TRUE)/sum(all_variants[WT == TRUE, ]$one_over_sigmasquared, 
            na.rm = TRUE)
        gr_lm <- dplyr::bind_cols(gr_lm, c(stop = stop1, wt = wt1))
        colnames(gr_lm)[i] <- paste0("block", i)
        all_variants$fitness_over_sigmasquared <- all_variants$fitness/(all_variants$sigma)^2
        all_variants$one_over_fitness_sigmasquared <- 1/(all_variants$sigma)^2
        stop1_fitness <- sum(all_variants[STOP == TRUE, ]$fitness_over_sigmasquared, na.rm = TRUE)/sum(all_variants[STOP == 
            TRUE, ]$one_over_fitness_sigmasquared, na.rm = TRUE)
        wt1_fitness <- sum(all_variants[WT == TRUE, ]$fitness_over_sigmasquared, na.rm = TRUE)/sum(all_variants[WT == TRUE, 
            ]$one_over_fitness_sigmasquared, na.rm = TRUE)
        fit_lm <- dplyr::bind_cols(fit_lm, c(stop = stop1_fitness, wt = wt1_fitness))
        colnames(fit_lm)[i] <- paste0("block", i)
        all_data[[i]] <- all_variants
    }
    all_data[[1]][, `:=`(nor_gr, growthrate)]
    all_data[[1]][, `:=`(nor_gr_sigma, growthrate_sigma)]
    all_data[[1]][, `:=`(nor_fitness, fitness)]
    all_data[[1]][, `:=`(nor_fitness_sigma, sigma)]
    for (i in 2:5) {
        formula1 <- as.formula(paste0("block1 ~ block", i))
        b1b2 <- summary(lm(formula = formula1, data = gr_lm))
        a2 <- b1b2$coefficients[2, 1]
        b2 <- b1b2$coefficients[1, 1]
        formula2 <- as.formula(paste0("block1 ~ block", i))
        c1c2 <- summary(lm(formula = formula2, data = fit_lm))
        d2 <- c1c2$coefficients[2, 1]
        e2 <- c1c2$coefficients[1, 1]
        all_data[[i]][, `:=`(nor_gr, growthrate * a2 + b2)]
        all_data[[i]][, `:=`(nor_gr_sigma, growthrate_sigma * a2)]
        all_data[[i]][, `:=`(nor_fitness, fitness * d2 + e2)]
        all_data[[i]][, `:=`(nor_fitness_sigma, sigma * d2)]
    }
    data_after_nor <- dplyr::bind_rows(all_data, .id = "block")
    return(data_after_nor)
}

