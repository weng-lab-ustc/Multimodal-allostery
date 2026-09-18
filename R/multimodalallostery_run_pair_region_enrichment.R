#' run pair region enrichment
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_run_pair_region_enrichment()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param pair_name Input or option used by this workflow; see Usage for the exact interface.
#' @param input_files Input or option used by this workflow; see Usage for the exact interface.
#' @param anno Input or option used by this workflow; see Usage for the exact interface.
#' @param structure_regions Input or option used by this workflow; see Usage for the exact interface.
#' @param fixed_threshold Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_sites_map Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_run_pair_region_enrichment <- function (pair_name, input_files, anno, structure_regions, fixed_threshold = 0.40000000000000002, binding_sites_map = list(RAF1 = c(21, 
    25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 
    73, 74), K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71), K13 = c(63, 68, 87, 88, 90, 91, 
    92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138), K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 
    97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137))) 
{
    assays <- strsplit(pair_name, " vs ")[[1]]
    x <- assays[1]
    y <- assays[2]
    cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
    cat("Analyzing:", pair_name, "\n")
    cat("Fixed threshold:", fixed_threshold, "kcal/mol\n")
    cat(paste0(rep("=", 60), collapse = ""), "\n")
    prepared <- multimodalallostery_prepare_region_merged_data_with_fdr(input_x = input_files[[x]], input_y = input_files[[y]], 
        assay_x = x, assay_y = y, anno = anno, fixed_threshold = fixed_threshold, binding_sites_map = binding_sites_map)
    df <- prepared$data
    threshold_x <- prepared$threshold_x
    threshold_y <- prepared$threshold_y
    df <- multimodalallostery_classify_two_step(merged_data = df, threshold_x = threshold_x, threshold_y = threshold_y, assay_x = x, 
        assay_y = y)
    df[, `:=`(category, data.table::fifelse(final_classification %in% c("Both promoting", "Both disrupting"), "Correlated", 
        data.table::fifelse(final_classification %in% c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y"), 
            "Anti-correlated", "Other")))]
    all_results <- list()
    for (region_name in names(structure_regions)) {
        region_residues <- structure_regions[[region_name]]
        plot_df <- df[, .(frac = mean(Pos_real %in% region_residues), n = .N, se = sqrt(mean(Pos_real %in% region_residues) * 
            (1 - mean(Pos_real %in% region_residues))/.N)), by = category]
        plot_df[, `:=`(region, region_name)]
        plot_df[, `:=`(pair, pair_name)]
        categories_to_test <- c("Correlated", "Anti-correlated", "Other")
        or_list <- lapply(categories_to_test, function(cat) {
            res <- multimodalallostery_calc_or_original(df = df, region_residues = region_residues, cat = cat)
            data.table::data.table(pair = pair_name, region = region_name, category = cat, OR = res$OR, OR_low = res$OR_low, 
                OR_high = res$OR_high, p = res$p)
        })
        or_df <- data.table::rbindlist(or_list)
        all_results[[region_name]] <- list(plot = plot_df, or = or_df)
    }
    combined_plot <- data.table::rbindlist(lapply(all_results, `[[`, "plot"))
    combined_or <- data.table::rbindlist(lapply(all_results, `[[`, "or"))
    cat("\n\u5206\u7c7b\u7edf\u8ba1:\n")
    cat("  Total mutations:", nrow(df), "\n")
    cat("  Correlated:", sum(df$category == "Correlated"), "\n")
    cat("  Anti-correlated:", sum(df$category == "Anti-correlated"), "\n")
    cat("  Other:", sum(df$category == "Other"), "\n")
    list(plot = combined_plot, or = combined_or, full_data = df)
}

