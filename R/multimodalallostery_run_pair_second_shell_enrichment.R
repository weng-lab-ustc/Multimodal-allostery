#' Run Pair Second Shell Enrichment.
#'
#' Reusable utility logic consolidated from the selected Figure/Panel implementations.
#'
#' @param pair_name Value supplied for `pair_name`.
#' @param input_files Value supplied for `input_files`.
#' @param anno Annotation data or annotation file path used by the analysis.
#' @param structure_regions Named list of structural regions and residue positions.
#' @param second_shell_map Named list mapping assays to second-shell residue positions.
#' @param binding_sites_map Named list mapping assays to binding-site residue positions.
#'
#' @return The transformed result produced by the function.
#' @export
multimodalallostery_run_pair_second_shell_enrichment <- function(pair_name, input_files, anno, structure_regions, second_shell_map, binding_sites_map) {
    assays <- strsplit(pair_name, " vs ")[[1]]
    x <- assays[1]
    y <- assays[2]
    cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
    cat("Analyzing:", pair_name, "\n")
    cat(paste0(rep("=", 60), collapse = ""), "\n")
    prepared <- multimodalallostery_prepare_region_merged_data_with_fdr(input_files[[x]], input_files[[y]], x, y, anno, binding_sites_map)
    df <- prepared$data
    threshold_x <- prepared$threshold_x
    threshold_y <- prepared$threshold_y
    df <- multimodalallostery_classify_two_step(df, threshold_x, threshold_y, x, y)
    df[, `:=`(category, fifelse(final_classification %in% c("Both promoting", "Both disrupting"), "Correlated", fifelse(final_classification %in% 
        c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y"), "Anti-correlated", fifelse(final_classification == 
        "Allosteric only in X", "Allosteric only in X", fifelse(final_classification == "Allosteric only in Y", "Allosteric only in Y", 
        "Other")))))]
    cat("\nClassification summary:\n")
    cat("  Total mutations:", nrow(df), "\n")
    cat("  Correlated:", sum(df$category == "Correlated"), "\n")
    cat("  Anti-correlated:", sum(df$category == "Anti-correlated"), "\n")
    cat("  Allosteric only in X:", sum(df$category == "Allosteric only in X"), "\n")
    cat("  Allosteric only in Y:", sum(df$category == "Allosteric only in Y"), "\n")
    cat("  Other:", sum(df$category == "Other"), "\n")
    all_results <- list()
    for (region_name in names(structure_regions)) {
        region_residues <- structure_regions[[region_name]]
        plot_df <- df[, .(frac = mean(Pos_real %in% region_residues), n = .N, se = sqrt(mean(Pos_real %in% region_residues) * 
            (1 - mean(Pos_real %in% region_residues))/.N)), by = category]
        plot_df[, `:=`(region, region_name)]
        plot_df[, `:=`(pair, pair_name)]
        categories_to_test <- c("Correlated", "Anti-correlated", "Allosteric only in X", "Allosteric only in Y", "Other")
        or_list <- lapply(categories_to_test, function(cat) {
            res <- multimodalallostery_calc_or_original(df, region_residues, cat)
            data.table(pair = pair_name, region = region_name, category = cat, OR = res$OR, OR_low = res$OR_low, OR_high = res$OR_high, 
                p = res$p)
        })
        or_df <- rbindlist(or_list)
        all_results[[region_name]] <- list(plot = plot_df, or = or_df)
    }
    cat("\n--- Second Shell Analysis ---\n")
    second_shell_x <- second_shell_map[[x]]
    second_shell_y <- second_shell_map[[y]]
    cat("  ", x, "second shell positions:", length(second_shell_x), "\n")
    cat("  ", y, "second shell positions:", length(second_shell_y), "\n")
    if (length(second_shell_x) > 0) {
        cat("\n  Analyzing", x, "second shell\n")
        plot_df_shell_x <- df[, .(frac = mean(Pos_real %in% second_shell_x), n = .N, se = sqrt(mean(Pos_real %in% second_shell_x) * 
            (1 - mean(Pos_real %in% second_shell_x))/.N)), by = category]
        plot_df_shell_x <- plot_df_shell_x[category %in% c("Allosteric only in X", "Other")]
        plot_df_shell_x[, `:=`(region, paste0("Second Shell (", x, ")"))]
        plot_df_shell_x[, `:=`(pair, pair_name)]
        categories_to_test <- c("Allosteric only in X", "Other")
        or_list_shell_x <- lapply(categories_to_test, function(cat) {
            res <- multimodalallostery_calc_or_original(df, second_shell_x, cat)
            data.table(pair = pair_name, region = paste0("Second Shell (", x, ")"), category = cat, OR = res$OR, OR_low = res$OR_low, 
                OR_high = res$OR_high, p = res$p)
        })
        or_df_shell_x <- rbindlist(or_list_shell_x)
        all_results[[paste0("Second Shell (", x, ")")]] <- list(plot = plot_df_shell_x, or = or_df_shell_x)
    }
    if (length(second_shell_y) > 0) {
        cat("\n  Analyzing", y, "second shell\n")
        plot_df_shell_y <- df[, .(frac = mean(Pos_real %in% second_shell_y), n = .N, se = sqrt(mean(Pos_real %in% second_shell_y) * 
            (1 - mean(Pos_real %in% second_shell_y))/.N)), by = category]
        plot_df_shell_y <- plot_df_shell_y[category %in% c("Allosteric only in Y", "Other")]
        plot_df_shell_y[, `:=`(region, paste0("Second Shell (", y, ")"))]
        plot_df_shell_y[, `:=`(pair, pair_name)]
        categories_to_test <- c("Allosteric only in Y", "Other")
        or_list_shell_y <- lapply(categories_to_test, function(cat) {
            res <- multimodalallostery_calc_or_original(df, second_shell_y, cat)
            data.table(pair = pair_name, region = paste0("Second Shell (", y, ")"), category = cat, OR = res$OR, OR_low = res$OR_low, 
                OR_high = res$OR_high, p = res$p)
        })
        or_df_shell_y <- rbindlist(or_list_shell_y)
        all_results[[paste0("Second Shell (", y, ")")]] <- list(plot = plot_df_shell_y, or = or_df_shell_y)
    }
    combined_plot <- rbindlist(lapply(all_results, `[[`, "plot"))
    combined_or <- rbindlist(lapply(all_results, `[[`, "or"))
    cat("\nFinal regions:", paste(unique(combined_plot$region), collapse = ", "), "\n")
    list(plot = combined_plot, or = combined_or, full_data = df)
}

