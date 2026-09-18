#' K13 19 Get Ob Pre Fitness Binding Correlation 3blocks.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param prediction Value supplied for `prediction`.
#' @param block1_dimsum_df Value supplied for `block1_dimsum_df`.
#' @param block2_dimsum_df Value supplied for `block2_dimsum_df`.
#' @param block3_dimsum_df Value supplied for `block3_dimsum_df`.
#' @param assay_sele Selected assay identifier.
#' @param wt_aa_input Value supplied for `wt_aa_input`.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_k13_19_get_ob_pre_fitness_binding_correlation_3blocks <- function(prediction = "path/to/predicted_phenotypes_all.txt", block1_dimsum_df = "path/to/block1_fitness.RData", 
    block2_dimsum_df = "path/to/block2_fitness.RData", block3_dimsum_df = "path/to/block3_fitness.RData", assay_sele = "RAF1", 
    wt_aa_input = wt_aa) {
    pre <- fread(prediction)
    pre_pos <- krasddpcams__pos_id(input = pre, wt_aa = wt_aa_input)
    load(block1_dimsum_df)
    block1 <- as.data.table(all_variants)
    load(block2_dimsum_df)
    block2 <- as.data.table(all_variants)
    load(block3_dimsum_df)
    block3 <- as.data.table(all_variants)
    data_before_nor <- rbind(block1 = block1, block2 = block2, block3 = block3, idcol = "block", fill = TRUE)
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
    multimodalallostery_calculate_wt_fitness <- function(block_name) {
        wt_fitness_block <- data_before_nor[WT == TRUE & block == block_name, ]
        wt_fitness <- sum(wt_fitness_block$fitness_over_sigmasquared, na.rm = TRUE)/sum(wt_fitness_block$one_over_fitness_sigmasquared, 
            na.rm = TRUE)
        return(wt_fitness)
    }
    wt1_fitness <- multimodalallostery_calculate_wt_fitness("block1")
    wt2_fitness <- multimodalallostery_calculate_wt_fitness("block2")
    wt3_fitness <- multimodalallostery_calculate_wt_fitness("block3")
    scaling_data_fitness <- data.frame(block1 = c(stop1_fitness, wt1_fitness), block2 = c(stop2_fitness, wt2_fitness), block3 = c(stop3_fitness, 
        wt3_fitness))
    multimodalallostery_calculate_scaling_params <- function(target_block) {
        lm_model <- lm(formula = block1 ~ get(target_block), data = scaling_data_fitness)
        return(list(slope = lm_model$coefficients[[2]], intercept = lm_model$coefficients[[1]]))
    }
    params_block2 <- multimodalallostery_calculate_scaling_params("block2")
    params_block3 <- multimodalallostery_calculate_scaling_params("block3")
    d2 <- params_block2$slope
    e2 <- params_block2$intercept
    d3 <- params_block3$slope
    e3 <- params_block3$intercept
    pre_nor <- pre_pos
    multimodalallostery_extract_prediction <- function(row) {
        return(row[78 + as.numeric(row[92])])
    }
    pre_nor$predicted_fitness <- apply(pre_nor, MARGIN = 1, FUN = multimodalallostery_extract_prediction)
    pre_nor$predicted_fitness <- as.numeric(pre_nor$predicted_fitness)
    assay_sele_df <- data.table(assay = c("K13", "K19", "K27", "K55", "PI3", "RAF1", "RAL", "SOS"), phenotype_base = c(6, 
        9, 12, 17, 22, 26, 29, 32))
    if (!assay_sele %in% assay_sele_df$assay) {
        stop("Assay '", assay_sele, "' not found. Available assays: ", paste(assay_sele_df$assay, collapse = ", "))
    }
    base_pheno <- assay_sele_df[assay == assay_sele, phenotype_base]
    pre_nor[phenotype == base_pheno, `:=`(pre_nor_mean_fitness, mean)]
    pre_nor[phenotype == base_pheno, `:=`(pre_nor_fitness_sigma, std)]
    pre_nor[phenotype == base_pheno, `:=`(ob_nor_fitness, fitness)]
    pre_nor[phenotype == base_pheno, `:=`(ob_nor_fitness_sigma, sigma)]
    pre_nor[phenotype == base_pheno, `:=`(pre_nor_fitness, predicted_fitness)]
    pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_mean_fitness, mean * d2 + e2)]
    pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_fitness_sigma, std * d2)]
    pre_nor[phenotype == base_pheno + 1, `:=`(ob_nor_fitness, fitness * d2 + e2)]
    pre_nor[phenotype == base_pheno + 1, `:=`(ob_nor_fitness_sigma, sigma * d2)]
    pre_nor[phenotype == base_pheno + 1, `:=`(pre_nor_fitness, predicted_fitness * d2 + e2)]
    pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_mean_fitness, mean * d3 + e3)]
    pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_fitness_sigma, std * d3)]
    pre_nor[phenotype == base_pheno + 2, `:=`(ob_nor_fitness, fitness * d3 + e3)]
    pre_nor[phenotype == base_pheno + 2, `:=`(ob_nor_fitness_sigma, sigma * d3)]
    pre_nor[phenotype == base_pheno + 2, `:=`(pre_nor_fitness, predicted_fitness * d3 + e3)]
    return(pre_nor)
}

