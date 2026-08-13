#' Run Pair Region Enrichment.
#'
#' Reusable utility logic consolidated from the selected Figure/Panel implementations.
#'
#' @param pair_name Value supplied for `pair_name`.
#' @param input_files Value supplied for `input_files`.
#' @param anno Annotation data or annotation file path used by the analysis.
#' @param structure_regions Named list of structural regions and residue positions.
#' @param binding_sites_map Named list mapping assays to binding-site residue positions.
#'
#' @return The transformed result produced by the function.
#' @export
multimodalallostery_run_pair_region_enrichment <- function(pair_name, input_files, anno, structure_regions, binding_sites_map) {
    assays <- strsplit(pair_name, " vs ")[[1]]
    x <- assays[1]
    y <- assays[2]
    prepared <- multimodalallostery_prepare_region_merged_data_with_fdr(input_files[[x]], input_files[[y]], x, y, anno, binding_sites_map)
    df <- prepared$data
    threshold_x <- prepared$threshold_x
    threshold_y <- prepared$threshold_y
    df <- multimodalallostery_classify_two_step(df, threshold_x, threshold_y, x, y)
    df[, `:=`(category, fifelse(final_classification %in% c("Both promoting", "Both disrupting"), "Correlated", fifelse(final_classification %in% 
        c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y"), "Anti-correlated", "Other")))]
    all_results <- list()
    for (region_name in names(structure_regions)) {
        region_residues <- structure_regions[[region_name]]
        plot_df <- df[, .(frac = mean(Pos_real %in% region_residues), n = .N, se = sqrt(mean(Pos_real %in% region_residues) * 
            (1 - mean(Pos_real %in% region_residues))/.N)), by = category]
        plot_df[, `:=`(region, region_name)]
        plot_df[, `:=`(pair, pair_name)]
        categories_to_test <- c("Correlated", "Anti-correlated", "Other")
        or_list <- lapply(categories_to_test, function(cat) {
            res <- multimodalallostery_calc_or_original(df, region_residues, cat)
            data.table(pair = pair_name, region = region_name, category = cat, OR = res$OR, OR_low = res$OR_low, OR_high = res$OR_high, 
                p = res$p)
        })
        or_df <- rbindlist(or_list)
        all_results[[region_name]] <- list(plot = plot_df, or = or_df)
    }
    combined_plot <- rbindlist(lapply(all_results, `[[`, "plot"))
    combined_or <- rbindlist(lapply(all_results, `[[`, "or"))
    cat("\n", pair_name, "\n")
    cat("  Total mutations:", nrow(df), "\n")
    cat("  Correlated:", sum(df$category == "Correlated"), "\n")
    cat("  Anti-correlated:", sum(df$category == "Anti-correlated"), "\n")
    cat("  Other:", sum(df$category == "Other"), "\n")
    list(plot = combined_plot, or = combined_or, full_data = df)
}

