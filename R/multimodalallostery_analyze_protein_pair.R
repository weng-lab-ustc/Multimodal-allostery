#' Analyze Protein Pair.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param protein_x Value supplied for `protein_x`.
#' @param protein_y Value supplied for `protein_y`.
#' @param input_file_x Value supplied for `input_file_x`.
#' @param input_file_y Value supplied for `input_file_y`.
#' @param anno_file Value supplied for `anno_file`.
#' @param binding_sites_map Named assay-to-binding-site residue mapping.
#' @param nbp_residues Residue positions assigned to the nucleotide-binding pocket.
#' @param legend_order Ordered classification labels used in downstream plots.
#' @param verbose Value supplied for `verbose`.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_analyze_protein_pair <- function(protein_x, protein_y, input_file_x, input_file_y, anno_file,
                                 binding_sites_map, nbp_residues, legend_order, verbose = TRUE) {
    prepared <- multimodalallostery_prepare_pairwise_merged_data_with_fdr(input_file_x, input_file_y, protein_x, protein_y, anno_file, binding_sites_map)
    merged_data <- prepared$data
    threshold_x <- prepared$threshold_x
    threshold_y <- prepared$threshold_y
    merged_data[, `:=`(pass_FDR_x, p_adj_x < 0.05)]
    merged_data[, `:=`(pass_FDR_y, p_adj_y < 0.05)]
    merged_data[, `:=`(direction_class, multimodalallostery_classify_by_direction(get(paste0("ddG_", protein_x)), get(paste0("ddG_", protein_y)), 
        threshold_x, threshold_y))]
    merged_data[, `:=`(final_classification, multimodalallostery_reclassify_by_fdr(direction_class, pass_FDR_x, pass_FDR_y))]
    merged_data[, `:=`(final_classification, factor(final_classification, levels = legend_order))]
    merged_data[, `:=`(is_NBP, FALSE)]
    merged_data[final_classification %in% c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y"), `:=`(is_NBP, 
        Pos_real %in% nbp_residues)]
    if (verbose) {
        cat("\n")
        cat(rep("=", 60), sep = "", collapse = "")
        cat("\n", protein_x, "vs", protein_y, "\u5206\u7C7B\u7ED3\u679C\n")
        cat(rep("=", 60), sep = "", collapse = "")
        cat("\n")
        cat("\nThreshold (", protein_x, "): ", threshold_x, " kcal/mol", sep = "")
        cat("\nThreshold (", protein_y, "): ", threshold_y, " kcal/mol", sep = "")
        cat("\nFDR Threshold: 0.05\n")
        both_promoting <- merged_data[final_classification == "Both promoting", mt]
        both_disrupting <- merged_data[final_classification == "Both disrupting", mt]
        promo_x_disrupt_y <- merged_data[final_classification == "Promoting in X / Disrupting in Y", mt]
        disrupt_x_promo_y <- merged_data[final_classification == "Disrupting in X / Promoting in Y", mt]
        only_x <- merged_data[final_classification == "Allosteric only in X", mt]
        only_y <- merged_data[final_classification == "Allosteric only in Y", mt]
        anticorrelated_nbp <- merged_data[final_classification %in% c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y") & 
            is_NBP == TRUE, .N]
        multimodalallostery_print_mutation_list(both_promoting, paste0("\n1. CORRELATED - Both promoting (", protein_x, " & ", protein_y, "):"))
        multimodalallostery_print_mutation_list(both_disrupting, paste0("\n2. CORRELATED - Both disrupting (", protein_x, " & ", protein_y, "):"))
        multimodalallostery_print_mutation_list(promo_x_disrupt_y, paste0("\n3. ANTICORRELATED - Promoting in ", protein_x, " / Disrupting in ", 
            protein_y, ":"))
        multimodalallostery_print_mutation_list(disrupt_x_promo_y, paste0("\n4. ANTICORRELATED - Disrupting in ", protein_x, " / Promoting in ", 
            protein_y, ":"))
        multimodalallostery_print_mutation_list(only_x, paste0("\n5. INDEPENDENT - Allosteric only in ", protein_x, ":"))
        multimodalallostery_print_mutation_list(only_y, paste0("\n6. INDEPENDENT - Allosteric only in ", protein_y, ":"))
        cat("\n", rep("=", 60), sep = "", collapse = "")
        cat("\nStatistical Summary:\n")
        cat("  Total number of mutations:", nrow(merged_data), "\n")
        cat("  Number of significant mutations:", sum(merged_data$final_classification != "Not significant (FDR >= 0.05)"), 
            "\n")
        cat("  Correlated (increase):", length(both_promoting), "\n")
        cat("  Correlated (inhibit):", length(both_disrupting), "\n")
        cat("  Anticorrelated (", protein_x, "increase/", protein_y, "inhibit):", length(promo_x_disrupt_y), "\n")
        cat("  Anticorrelated (", protein_x, "inhibit/", protein_y, "increase):", length(disrupt_x_promo_y), "\n")
        cat("    - Among them, the NBP site:", anticorrelated_nbp, "\n")
        cat("  Independent (only", protein_x, "):", length(only_x), "\n")
        cat("  Independent (only", protein_y, "):", length(only_y), "\n")
        cat("  Not significant (FDR >= 0.05):", sum(merged_data$final_classification == "Not significant (FDR >= 0.05)"), 
            "\n")
        cat(rep("=", 60), sep = "", collapse = "")
        cat("\n")
        cat("\nComparative Classification (Direction vs. Outcome):\n")
        comparison <- merged_data[, .N, by = .(direction_class, final_classification)]
        print(comparison)
        mutations_list <- list(both_promoting = both_promoting, both_disrupting = both_disrupting, promo_x_disrupt_y = promo_x_disrupt_y, 
            disrupt_x_promo_y = disrupt_x_promo_y, only_x = only_x, only_y = only_y)
    }
    else {
        mutations_list <- NULL
    }
    return(list(data = merged_data, thresholds = c(threshold_x, threshold_y), names = c(protein_x, protein_y), mutations = mutations_list))
}

