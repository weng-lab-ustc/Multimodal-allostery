#' Merge Ddgf Fitness Blocks.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param prediction Value supplied for `prediction`.
#' @param folding_ddG Value supplied for `folding_ddG`.
#' @param block1_dimsum_df Value supplied for `block1_dimsum_df`.
#' @param block2_dimsum_df Value supplied for `block2_dimsum_df`.
#' @param block3_dimsum_df Value supplied for `block3_dimsum_df`.
#' @param block4_dimsum_df Value supplied for `block4_dimsum_df`.
#' @param block5_dimsum_df Value supplied for `block5_dimsum_df`.
#' @param wt_aa_input Value supplied for `wt_aa_input`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_merge_ddgf_fitness_blocks <- function(prediction = prediction, folding_ddG = folding_ddG, block1_dimsum_df = block1_dimsum_df, 
    block2_dimsum_df = block2_dimsum_df, block3_dimsum_df = block3_dimsum_df, block4_dimsum_df = block4_dimsum_df, block5_dimsum_df = block5_dimsum_df, 
    wt_aa_input = wt_aa_input) {
    pre <- fread(prediction)
    folding_ddG <- fread(folding_ddG)
    pre_pos <- krasddpcams__pos_id(input = pre, wt_aa = wt_aa_input)
    load(block1_dimsum_df)
    block1 <- as.data.table(all_variants)
    load(block2_dimsum_df)
    block2 <- as.data.table(all_variants)
    load(block3_dimsum_df)
    block3 <- as.data.table(all_variants)
    load(block4_dimsum_df)
    block4 <- as.data.table(all_variants)
    load(block5_dimsum_df)
    block5 <- as.data.table(all_variants)
    data_before_nor <- rbind(block1 = block1, block2 = block2, block3 = block3, block4 = block4, block5 = block5, idcol = "block", 
        fill = TRUE)
    data_before_nor$fitness_over_sigmasquared <- data_before_nor$fitness/(data_before_nor$sigma)^2
    data_before_nor$one_over_fitness_sigmasquared <- 1/(data_before_nor$sigma)^2
    multimodalallostery_calculate_stop_fitness <- function(block_name) {
        dead_fitness <- data_before_nor[STOP == TRUE & block == block_name, ]
        stop_fitness <- sum(dead_fitness$fitness_over_sigmasquared, na.rm = TRUE)/sum(dead_fitness$one_over_fitness_sigmasquared, 
            na.rm = TRUE)
        return(stop_fitness)
    }
    stop1_fitness <- multimodalallostery_calculate_stop_fitness("block1")
    stop2_fitness <- multimodalallostery_calculate_stop_fitness("block2")
    stop3_fitness <- multimodalallostery_calculate_stop_fitness("block3")
    stop4_fitness <- multimodalallostery_calculate_stop_fitness("block4")
    stop5_fitness <- multimodalallostery_calculate_stop_fitness("block5")
    multimodalallostery_calculate_wt_fitness <- function(block_name) {
        wt_fitness_block <- data_before_nor[WT == TRUE & block == block_name, ]
        wt_fitness <- sum(wt_fitness_block$fitness_over_sigmasquared, na.rm = TRUE)/sum(wt_fitness_block$one_over_fitness_sigmasquared, 
            na.rm = TRUE)
        return(wt_fitness)
    }
    wt1_fitness <- multimodalallostery_calculate_wt_fitness("block1")
    wt2_fitness <- multimodalallostery_calculate_wt_fitness("block2")
    wt3_fitness <- multimodalallostery_calculate_wt_fitness("block3")
    wt4_fitness <- multimodalallostery_calculate_wt_fitness("block4")
    wt5_fitness <- multimodalallostery_calculate_wt_fitness("block5")
    scaling_data_fitness <- data.frame(block1 = c(stop1_fitness, wt1_fitness), block2 = c(stop2_fitness, wt2_fitness), block3 = c(stop3_fitness, 
        wt3_fitness), block4 = c(stop4_fitness, wt4_fitness), block5 = c(stop5_fitness, wt5_fitness))
    multimodalallostery_calculate_scaling_params <- function(target_block) {
        lm_model <- lm(formula = block1 ~ get(target_block), data = scaling_data_fitness)
        return(list(slope = lm_model$coefficients[[2]], intercept = lm_model$coefficients[[1]]))
    }
    params_block2 <- multimodalallostery_calculate_scaling_params("block2")
    params_block3 <- multimodalallostery_calculate_scaling_params("block3")
    params_block4 <- multimodalallostery_calculate_scaling_params("block4")
    params_block5 <- multimodalallostery_calculate_scaling_params("block5")
    d2 <- params_block2$slope
    e2 <- params_block2$intercept
    d3 <- params_block3$slope
    e3 <- params_block3$intercept
    d4 <- params_block4$slope
    e4 <- params_block4$intercept
    d5 <- params_block5$slope
    e5 <- params_block5$intercept
    pre_nor <- pre_pos
    multimodalallostery_extract_prediction <- function(row) {
        return(row[78 + as.numeric(row[92])])
    }
    pre_nor$predicted_fitness <- apply(pre_nor, MARGIN = 1, FUN = multimodalallostery_extract_prediction)
    pre_nor$predicted_fitness <- as.numeric(pre_nor$predicted_fitness)
    multimodalallostery_extract_additive_trait0 <- function(row) {
        return(row[92 + as.numeric(row[92]) * 2 - 1])
    }
    multimodalallostery_extract_additive_trait1 <- function(row) {
        return(row[92 + as.numeric(row[92]) * 2])
    }
    pre_nor$additive_trait0 <- apply(pre_nor, MARGIN = 1, FUN = multimodalallostery_extract_additive_trait0)
    pre_nor$additive_trait0 <- as.numeric(pre_nor$additive_trait0)
    pre_nor$additive_trait1 <- apply(pre_nor, MARGIN = 1, FUN = multimodalallostery_extract_additive_trait1)
    pre_nor$additive_trait1 <- as.numeric(pre_nor$additive_trait1)
    pre_nor[, `:=`(additive_trait, additive_trait0 + additive_trait1)]
    pre_nor[phenotype == 1, `:=`(pre_nor_mean_fitness, mean)]
    pre_nor[phenotype == 2, `:=`(pre_nor_mean_fitness, mean * d2 + e2)]
    pre_nor[phenotype == 3, `:=`(pre_nor_mean_fitness, mean * d3 + e3)]
    pre_nor[phenotype == 4, `:=`(pre_nor_mean_fitness, mean * d4 + e4)]
    pre_nor[phenotype == 5, `:=`(pre_nor_mean_fitness, mean * d5 + e5)]
    pre_nor[phenotype == 1, `:=`(pre_nor_fitness_sigma, std)]
    pre_nor[phenotype == 2, `:=`(pre_nor_fitness_sigma, std * d2)]
    pre_nor[phenotype == 3, `:=`(pre_nor_fitness_sigma, std * d3)]
    pre_nor[phenotype == 4, `:=`(pre_nor_fitness_sigma, std * d4)]
    pre_nor[phenotype == 5, `:=`(pre_nor_fitness_sigma, std * d5)]
    pre_nor[phenotype == 1, `:=`(ob_nor_fitness, fitness)]
    pre_nor[phenotype == 2, `:=`(ob_nor_fitness, fitness * d2 + e2)]
    pre_nor[phenotype == 3, `:=`(ob_nor_fitness, fitness * d3 + e3)]
    pre_nor[phenotype == 4, `:=`(ob_nor_fitness, fitness * d4 + e4)]
    pre_nor[phenotype == 5, `:=`(ob_nor_fitness, fitness * d5 + e5)]
    pre_nor[phenotype == 1, `:=`(ob_nor_fitness_sigma, sigma)]
    pre_nor[phenotype == 2, `:=`(ob_nor_fitness_sigma, sigma * d2)]
    pre_nor[phenotype == 3, `:=`(ob_nor_fitness_sigma, sigma * d3)]
    pre_nor[phenotype == 4, `:=`(ob_nor_fitness_sigma, sigma * d4)]
    pre_nor[phenotype == 5, `:=`(ob_nor_fitness_sigma, sigma * d5)]
    pre_nor[phenotype == 1, `:=`(pre_nor_fitness, predicted_fitness)]
    pre_nor[phenotype == 2, `:=`(pre_nor_fitness, predicted_fitness * d2 + e2)]
    pre_nor[phenotype == 3, `:=`(pre_nor_fitness, predicted_fitness * d3 + e3)]
    pre_nor[phenotype == 4, `:=`(pre_nor_fitness, predicted_fitness * d4 + e4)]
    pre_nor[phenotype == 5, `:=`(pre_nor_fitness, predicted_fitness * d5 + e5)]
    return(pre_nor)
}

