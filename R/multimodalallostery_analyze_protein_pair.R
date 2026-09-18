#' analyze protein pair
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_analyze_protein_pair()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param protein_x Input or option used by this workflow; see Usage for the exact interface.
#' @param protein_y Input or option used by this workflow; see Usage for the exact interface.
#' @param input_file_x Input or option used by this workflow; see Usage for the exact interface.
#' @param input_file_y Input or option used by this workflow; see Usage for the exact interface.
#' @param anno_file Input or option used by this workflow; see Usage for the exact interface.
#' @param fixed_threshold Input or option used by this workflow; see Usage for the exact interface.
#' @param verbose Input or option used by this workflow; see Usage for the exact interface.
#' @param legend_order Input or option used by this workflow; see Usage for the exact interface.
#' @param nbp_residues Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_sites_map Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_analyze_protein_pair <- function (protein_x, protein_y, input_file_x, input_file_y, anno_file, fixed_threshold = 0.40000000000000002, verbose = TRUE, 
    legend_order = c("Both promoting", "Both disrupting", "Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y", 
        "Allosteric only in X", "Allosteric only in Y", "Other (neutral in both)", "Not significant (FDR >= 0.05)"), nbp_residues = c(12, 
        13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147), binding_sites_map = list(RAF1 = c(21, 
        25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 
        70, 73, 74), K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71), K13 = c(63, 68, 87, 88, 
        90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138), K19 = c(68, 87, 88, 90, 91, 
        92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137))) 
{
    prepared <- multimodalallostery_prepare_pairwise_merged_data_with_fdr(input_x = input_file_x, input_y = input_file_y, 
        assay_x = protein_x, assay_y = protein_y, anno = anno_file, fixed_threshold = fixed_threshold, binding_sites_map = binding_sites_map)
    merged_data <- prepared$data
    threshold_x <- prepared$threshold_x
    threshold_y <- prepared$threshold_y
    merged_data[, `:=`(pass_FDR_x, p_adj_x < 0.050000000000000003)]
    merged_data[, `:=`(pass_FDR_y, p_adj_y < 0.050000000000000003)]
    merged_data[, `:=`(direction_class, multimodalallostery_classify_by_direction(ddG_x = get(paste0("ddG_", protein_x)), 
        ddG_y = get(paste0("ddG_", protein_y)), threshold_x = threshold_x, threshold_y = threshold_y))]
    merged_data[, `:=`(final_classification, multimodalallostery_reclassify_by_fdr(direction_class = direction_class, pass_FDR_x = pass_FDR_x, 
        pass_FDR_y = pass_FDR_y))]
    merged_data[, `:=`(final_classification, factor(final_classification, levels = legend_order))]
    merged_data[, `:=`(is_NBP, FALSE)]
    merged_data[final_classification %in% c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y"), `:=`(is_NBP, 
        Pos_real %in% nbp_residues)]
    if (verbose) {
        cat("\n")
        cat(rep("=", 60), sep = "", collapse = "")
        cat("\n", protein_x, "vs", protein_y, "\u5206\u7c7b\u7ed3\u679c\n")
        cat(rep("=", 60), sep = "", collapse = "")
        cat("\n")
        cat("\n\u56fa\u5b9a\u9608\u503c: ", fixed_threshold, " kcal/mol (\u7528\u4e8e", protein_x, "\u548c", protein_y, ")", sep = "")
        cat("\nFDR\u9608\u503c: 0.05\n")
        both_promoting <- merged_data[final_classification == "Both promoting", mt]
        both_disrupting <- merged_data[final_classification == "Both disrupting", mt]
        promo_x_disrupt_y <- merged_data[final_classification == "Promoting in X / Disrupting in Y", mt]
        disrupt_x_promo_y <- merged_data[final_classification == "Disrupting in X / Promoting in Y", mt]
        only_x <- merged_data[final_classification == "Allosteric only in X", mt]
        only_y <- merged_data[final_classification == "Allosteric only in Y", mt]
        anticorrelated_nbp <- merged_data[final_classification %in% c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y") & 
            is_NBP == TRUE, .N]
        multimodalallostery_print_mutation_list(mutations = both_promoting, title = paste0("\n1. CORRELATED - Both promoting (", 
            protein_x, " & ", protein_y, "):"))
        multimodalallostery_print_mutation_list(mutations = both_disrupting, title = paste0("\n2. CORRELATED - Both disrupting (", 
            protein_x, " & ", protein_y, "):"))
        multimodalallostery_print_mutation_list(mutations = promo_x_disrupt_y, title = paste0("\n3. ANTICORRELATED - Promoting in ", 
            protein_x, " / Disrupting in ", protein_y, ":"))
        multimodalallostery_print_mutation_list(mutations = disrupt_x_promo_y, title = paste0("\n4. ANTICORRELATED - Disrupting in ", 
            protein_x, " / Promoting in ", protein_y, ":"))
        multimodalallostery_print_mutation_list(mutations = only_x, title = paste0("\n5. INDEPENDENT - Allosteric only in ", 
            protein_x, ":"))
        multimodalallostery_print_mutation_list(mutations = only_y, title = paste0("\n6. INDEPENDENT - Allosteric only in ", 
            protein_y, ":"))
        cat("\n", rep("=", 60), sep = "", collapse = "")
        cat("\n\u7edf\u8ba1\u6c47\u603b:\n")
        cat("  \u603b\u7a81\u53d8\u6570:", nrow(merged_data), "\n")
        cat("  \u663e\u8457\u7a81\u53d8\u6570:", sum(merged_data$final_classification != "Not significant (FDR >= 0.05)"), "\n")
        cat("  Correlated (\u4fc3\u8fdb):", length(both_promoting), "\n")
        cat("  Correlated (\u7834\u574f):", length(both_disrupting), "\n")
        cat("  Anticorrelated (", protein_x, "\u4fc3\u8fdb/", protein_y, "\u7834\u574f):", length(promo_x_disrupt_y), "\n")
        cat("  Anticorrelated (", protein_x, "\u7834\u574f/", protein_y, "\u4fc3\u8fdb):", length(disrupt_x_promo_y), "\n")
        cat("    - \u5176\u4e2dNBP\u4f4d\u70b9:", anticorrelated_nbp, "\n")
        cat("  Independent (\u4ec5", protein_x, "):", length(only_x), "\n")
        cat("  Independent (\u4ec5", protein_y, "):", length(only_y), "\n")
        cat("  Not significant (FDR >= 0.05):", sum(merged_data$final_classification == "Not significant (FDR >= 0.05)"), 
            "\n")
        cat(rep("=", 60), sep = "", collapse = "")
        cat("\n")
        cat("\n\u5206\u7c7b\u5bf9\u6bd4\uff08\u65b9\u5411 vs \u6700\u7ec8\uff09:\n")
        comparison <- merged_data[, .N, by = .(direction_class, final_classification)]
        print(comparison)
        mutations_list <- list(both_promoting = both_promoting, both_disrupting = both_disrupting, promo_x_disrupt_y = promo_x_disrupt_y, 
            disrupt_x_promo_y = disrupt_x_promo_y, only_x = only_x, only_y = only_y)
    }
    else {
        mutations_list <- NULL
    }
    return(list(data = merged_data, thresholds = c(threshold_x, threshold_y), names = c(protein_x, protein_y), mutations = mutations_list, 
        fixed_threshold = fixed_threshold))
}

