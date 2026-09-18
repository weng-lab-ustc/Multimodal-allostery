#' Run Pair Second Shell Enrichment (Step-by-step 版本一致)
#'
#' @param pair_name Value supplied for `pair_name`.
#' @param input_files Value supplied for `input_files`.
#' @param anno Annotation data or annotation file path used by the analysis.
#' @param structure_regions Named list of structural regions and residue positions.
#' @param second_shell_map Named list mapping assays to second-shell residue positions.
#' @param binding_sites_map Named list mapping assays to binding-site residue positions.
#' @param fixed_threshold Fixed ddG threshold used for significance calls.
#'
#' @return The transformed result produced by the function.
#' @export
multimodalallostery_run_pair_second_shell_enrichment <- function(pair_name, input_files, anno,
                                                                structure_regions, second_shell_map,
                                                                binding_sites_map, fixed_threshold = 0.4) {
  assays <- strsplit(pair_name, " vs ")[[1]]; x <- assays[1]; y <- assays[2]
  cat("\n", strrep("=", 60), "\nAnalyzing:", pair_name, "\n", strrep("=", 60), "\n")
  prepared <- multimodalallostery_prepare_region_merged_data_with_fdr(input_x = input_files[[x]], input_y = input_files[[y]], assay_x = x, assay_y = y,
                                                                      anno = anno, fixed_threshold = fixed_threshold, binding_sites_map = binding_sites_map)
  df <- prepared$data
  df <- multimodalallostery_classify_two_step(df, prepared$threshold_x, prepared$threshold_y, x, y)
  df[, category := fifelse(final_classification %in% c("Both promoting", "Both disrupting"), "Correlated",
                    fifelse(final_classification %in% c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y"), "Anti-correlated",
                    fifelse(final_classification == "Allosteric only in X", "Allosteric only in X",
                    fifelse(final_classification == "Allosteric only in Y", "Allosteric only in Y", "Other"))))]
  cat("\nClassification summary:\n")
  for (k in c("Total mutations" = nrow(df), "Correlated" = sum(df$category == "Correlated"),
              "Anti-correlated" = sum(df$category == "Anti-correlated"),
              "Allosteric only in X" = sum(df$category == "Allosteric only in X"),
              "Allosteric only in Y" = sum(df$category == "Allosteric only in Y"),
              "Other" = sum(df$category == "Other")))
  cat("  ", names(k), ":", k, "\n")
  build_region <- function(region_name, region_residues, cats_plot, cats_or) {
    plot_df <- df[, .(frac = mean(Pos_real %in% region_residues), n = .N,
                      se = sqrt(mean(Pos_real %in% region_residues) *
                                (1 - mean(Pos_real %in% region_residues)) / .N)), by = category]
    plot_df <- plot_df[category %in% cats_plot]
    plot_df[, `:=`(region = region_name, pair = pair_name)]
    or_df <- rbindlist(lapply(cats_or, function(cat) {
      res <- multimodalallostery_calc_or_original(df, region_residues, cat)
      data.table(pair = pair_name, region = region_name, category = cat,
                 OR = res$OR, OR_low = res$OR_low, OR_high = res$OR_high, p = res$p)}))
    list(plot = plot_df, or = or_df) }
  all_cats <- c("Correlated", "Anti-correlated", "Allosteric only in X", "Allosteric only in Y", "Other")
  all_results <- lapply(names(structure_regions), function(rn)
    build_region(rn, structure_regions[[rn]], all_cats, all_cats))
  names(all_results) <- names(structure_regions)
  cat("\n--- Second Shell Analysis ---\n")
  second_shell_x <- second_shell_map[[x]]; second_shell_y <- second_shell_map[[y]]
  cat("  ", x, "second shell positions:", length(second_shell_x), "\n")
  cat("  ", y, "second shell positions:", length(second_shell_y), "\n")
  if (length(second_shell_x) > 0) { cat("\n  Analyzing", x, "second shell\n")
    all_results[[paste0("Second Shell (", x, ")")]] <-
      build_region(paste0("Second Shell (", x, ")"), second_shell_x,
                   c("Allosteric only in X", "Other"), c("Allosteric only in X", "Other")) }
  if (length(second_shell_y) > 0) { cat("\n  Analyzing", y, "second shell\n")
    all_results[[paste0("Second Shell (", y, ")")]] <-
      build_region(paste0("Second Shell (", y, ")"), second_shell_y,
                   c("Allosteric only in Y", "Other"), c("Allosteric only in Y", "Other")) }
  combined_plot <- rbindlist(lapply(all_results, `[[`, "plot"), fill = TRUE)
  combined_or   <- rbindlist(lapply(all_results, `[[`, "or"),   fill = TRUE)
  cat("\nFinal regions:", paste(unique(combined_plot$region), collapse = ", "), "\n")
  categories_to_plot <- c("Allosteric only in X", "Allosteric only in Y", "Other")
  color_map <- c("Allosteric only in X" = "#C68EFD", "Allosteric only in Y" = "#09B636", "Other" = "grey80")
  legend_labels <- c(paste0("Allosteric only in ", x), paste0("Allosteric only in ", y), "Other")
  plot_subset <- copy(combined_plot)[category %in% categories_to_plot]
  region_levels <- c(names(structure_regions), sort(grep("Second Shell", unique(plot_subset$region), value = TRUE)))
  region_levels <- region_levels[region_levels %in% unique(plot_subset$region)]
  plot_subset[, `:=`(region = factor(region, levels = region_levels), category = factor(category, levels = categories_to_plot))]
  or_subset <- copy(combined_or)[pair == pair_name & category %in% categories_to_plot & region %in% region_levels]
  or_subset[, `:=`(region = factor(region, levels = region_levels), x_offset = dplyr::case_when(category == "Allosteric only in X" ~ -0.25,
                                               category == "Allosteric only in Y" ~  0,
                                               category == "Other"               ~  0.25))]
  fixed_y_position <- 1.2
  max_y <- max(max(plot_subset$frac + plot_subset$se, na.rm = TRUE) * 1.15, fixed_y_position + 0.1, 1.0)
  ##or_subset[, label := paste0("OR = ", round(OR, 2), ifelse(p < 0.05, "*", ""))]
  or_subset[, label := paste0("OR = ", round(OR, 2))]
  p <- ggplot(plot_subset, aes(region, frac, fill = category)) +
    geom_col(position = position_dodge(0.9), width = 0.7) +
    #geom_errorbar(aes(ymin = frac - se, ymax = frac + se),  position = position_dodge(0.9), width = 0.15, size = 0.8, color = "#F1DD10", alpha = 0.8, na.rm = TRUE) +
    #geom_errorbar(aes(ymin = frac - se, ymax = frac - se), position = position_dodge(0.9), width = 0.3, size = 0.8, color = "#F1DD10", alpha = 0.8, na.rm = TRUE) +
    #geom_errorbar(aes(ymin = frac + se, ymax = frac + se), position = position_dodge(0.9), width = 0.3, size = 0.8, color = "#F1DD10", alpha = 0.8, na.rm = TRUE) +
    geom_text(data = or_subset, aes(x = as.numeric(region) + x_offset, y = fixed_y_position, label = label, color = category), size = 3.5, angle = 45, hjust = 0.5, vjust = 0) +
    scale_fill_manual(values = color_map, labels = legend_labels, drop = FALSE) +
    scale_color_manual(values = color_map, labels = legend_labels, guide = "none") +
    scale_y_continuous(limits = c(0, max_y), expand = expansion(mult = c(0, 0.02)), breaks = seq(0, 1, 0.2)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.5, size = 0.8) +
    labs(y = "Fraction of mutations in region", x = "Structural region") +
    theme_classic(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 15),
          axis.line = element_line(color = "black", size = 0.5),
          axis.ticks = element_line(color = "black", size = 0.5),
          panel.grid = element_blank(), legend.title = element_blank(),
          legend.position = "bottom",
          plot.margin = margin(t = 40, r = 10, b = 10, l = 10))
  list(plot = p, or = combined_or, fractions = combined_plot, full_data = df)
}