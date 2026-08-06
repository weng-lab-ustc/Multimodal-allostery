library(data.table)
library(ggplot2)
library(dplyr)

# =========================================================
PAIR_ORDER <- c(
  "RAF1 vs K55",
  "RAF1 vs K27",
  "K55 vs K27",
  "RAF1 vs K13",
  "RAF1 vs K19",
  "K55 vs K13",
  "K55 vs K19",
  "K27 vs K13",
  "K27 vs K19",
  "K13 vs K19"
)

# 
binding_sites_map <- list(
  RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
  K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74),
  K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71),
  K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138),
  K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
)

# =========================================================
core_residues <- c(4,6,7,8,9,10,11,14,15,16,17,18,19,20,21,22,23,24,40,
                   42,44,46,51,52,53,54,55,56,57,58,68,71,72,75,77,78,79,
                   80,81,82,83,84,89,90,92,93,96,97,99,100,101,103,109,110,
                   111,112,113,114,115,116,118,125,130,133,134,137,139,141,
                   142,143,144,145,146,151,152,155,156,157,158,159,160,162,163)

NBP <- c(12,13,14,15,16,17,18,28,29,30,32,34,35,57,60,61,116,117,119,120,145,146,147)

functional_loop_residues <- c(
  10, 11, 12, 13, 14, 15, 16, 17, 
  25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 
  58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76
)

beta_sheets <- c(3:9, 38:44, 51:57, 77:84, 109:115, 139:143)

# 
surface_residues <- c(1,2,3,5,12,13,25:39,41,43,45,47:50,59:67,
                      69,70,73,74,76,85:88,91,94,95,98,102,104:108,
                      117,119:124,126:129,131,132,135,136,138,140,
                      147:150,153,154,161,164:188)

alpha_helices <- c(15:24, 67:73, 87:104, 127:136, 148:166)

# =========================================================

load_mutation_data <- function(input, assay_sele) {
  if (is.character(input)) {
    ddG <- fread(input)
  } else {
    ddG <- as.data.table(input)
  }
  
  ddG[, Pos_real := Pos_ref + 1]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, mt := paste0(wt_codon, Pos_real, mt_codon)]
  ddG <- ddG[id != "WT"]
  
  result <- ddG[, .(mt, Pos_real, wt_codon, mt_codon, 
                    `mean_kcal/mol`, `std_kcal/mol`)]
  setnames(result, "mean_kcal/mol", "ddG")
  setnames(result, "std_kcal/mol", "ddG_std")
  result[, assay := assay_sele]
  
  return(result)
}

# =========================================================
classify_by_direction <- function(ddG_x, ddG_y, threshold_x, threshold_y) {
  sig_x <- abs(ddG_x) > threshold_x
  sig_y <- abs(ddG_y) > threshold_y
  
  disrupt_x <- ddG_x > threshold_x
  promote_x <- ddG_x < -threshold_x
  disrupt_y <- ddG_y > threshold_y
  promote_y <- ddG_y < -threshold_y
  
  result <- rep("neutral", length(ddG_x))
  
  #
  result[sig_x & sig_y & promote_x & promote_y] <- "Both promoting"
  result[sig_x & sig_y & disrupt_x & disrupt_y] <- "Both disrupting"
  result[sig_x & sig_y & promote_x & disrupt_y] <- "Promoting in X / Disrupting in Y"
  result[sig_x & sig_y & disrupt_x & promote_y] <- "Disrupting in X / Promoting in Y"
  
  #
  result[sig_x & !sig_y] <- "Allosteric only in X"
  result[!sig_x & sig_y] <- "Allosteric only in Y"
  
  #
  result[!sig_x & !sig_y] <- "Not significant (FDR >= 0.05)"
  
  return(result)
}

# =========================================================

reclassify_by_FDR <- function(direction_class, pass_FDR_x, pass_FDR_y) {
  result <- direction_class
  
  for (i in 1:length(direction_class)) {
    
    # 
    if (direction_class[i] %in% c("Both promoting", "Both disrupting", 
                                  "Promoting in X / Disrupting in Y", 
                                  "Disrupting in X / Promoting in Y")) {
      if (!(pass_FDR_x[i] & pass_FDR_y[i])) {
        # 
        if (pass_FDR_x[i] & !pass_FDR_y[i]) {
          result[i] <- "Allosteric only in X"
        } else if (!pass_FDR_x[i] & pass_FDR_y[i]) {
          result[i] <- "Allosteric only in Y"
        } else {
          result[i] <- "Not significant (FDR >= 0.05)"
        }
      }
    }
    
    # 
    else if (direction_class[i] == "Allosteric only in X") {
      if (!pass_FDR_x[i]) {
        result[i] <- "Not significant (FDR >= 0.05)"
      }
    }
    
    # 
    else if (direction_class[i] == "Allosteric only in Y") {
      if (!pass_FDR_y[i]) {
        result[i] <- "Not significant (FDR >= 0.05)"
      }
    }
  }
  
  return(result)
}

# =========================================================
calculate_threshold <- function(data, assay_sele, anno) {
  anno_data <- fread(anno)
  
  site_ddG <- data[, .(
    mean_ddG = sum(abs(ddG) / ddG_std^2, na.rm = TRUE) / sum(1 / ddG_std^2, na.rm = TRUE)
  ), by = Pos_real]
  
  site_sigma <- data[, .(
    sigma = sqrt(1 / sum(1 / ddG_std^2, na.rm = TRUE))
  ), by = Pos_real]
  
  site_stats <- merge(site_ddG, site_sigma, by = "Pos_real")
  site_stats <- merge(site_stats, anno_data, by.x = "Pos_real", by.y = "Pos", all.x = TRUE)
  
  scHAmin_col <- paste0("scHAmin_ligand_", assay_sele)
  binding_sites <- site_stats[get(scHAmin_col) < 5, Pos_real]
  
  threshold <- site_stats[Pos_real %in% binding_sites, 
                          sum(abs(mean_ddG) / sigma^2, na.rm = TRUE) / sum(1 / sigma^2, na.rm = TRUE)]
  
  return(threshold)
}

# =========================================================
prepare_merged_data_with_FDR <- function(input_x, input_y, assay_x, assay_y, anno) {
  data_x <- load_mutation_data(input_x, assay_x)
  data_y <- load_mutation_data(input_y, assay_y)
  
  threshold_x <- calculate_threshold(data_x, assay_x, anno)
  threshold_y <- calculate_threshold(data_y, assay_y, anno)
  
  if (!is.null(binding_sites_map[[assay_x]])) {
    data_x <- data_x[!(Pos_real %in% binding_sites_map[[assay_x]])]
  }
  if (!is.null(binding_sites_map[[assay_y]])) {
    data_y <- data_y[!(Pos_real %in% binding_sites_map[[assay_y]])]
  }
  
  data_x_clean <- data_x[, .(mt, Pos_real, ddG, ddG_std)]
  setnames(data_x_clean, "ddG", paste0("ddG_", assay_x))
  setnames(data_x_clean, "ddG_std", paste0("std_", assay_x))
  
  data_y_clean <- data_y[, .(mt, Pos_real, ddG, ddG_std)]
  setnames(data_y_clean, "ddG", paste0("ddG_", assay_y))
  setnames(data_y_clean, "ddG_std", paste0("std_", assay_y))
  
  merged_data <- merge(data_x_clean, data_y_clean, by = c("mt", "Pos_real"))
  
  pvalue_threshold <- function(av, se, threshold){
    zscore <- (abs(av) - threshold) / se
    2 * pnorm(abs(zscore), lower.tail = FALSE)
  }
  
  merged_data[, p_x := pvalue_threshold(get(paste0("ddG_", assay_x)),
                                        get(paste0("std_", assay_x)),
                                        threshold_x)]
  
  merged_data[, p_y := pvalue_threshold(get(paste0("ddG_", assay_y)),
                                        get(paste0("std_", assay_y)),
                                        threshold_y)]
  
  merged_data[, p_adj_x := p.adjust(p_x, "BH")]
  merged_data[, p_adj_y := p.adjust(p_y, "BH")]
  
  list(
    data = merged_data,
    threshold_x = threshold_x,
    threshold_y = threshold_y
  )
}

# =========================================================

classify_two_step <- function(merged_data, threshold_x, threshold_y, assay_x, assay_y) {
  
  ddG_x_col <- paste0("ddG_", assay_x)
  ddG_y_col <- paste0("ddG_", assay_y)
  
  #
  merged_data[, direction_class := classify_by_direction(
    get(ddG_x_col), get(ddG_y_col), threshold_x, threshold_y
  )]
  
  #
  merged_data[, pass_FDR_x := p_adj_x < 0.05]
  merged_data[, pass_FDR_y := p_adj_y < 0.05]
  
  #
  merged_data[, final_classification := reclassify_by_FDR(direction_class, pass_FDR_x, pass_FDR_y)]
  
  return(merged_data)
}

# =========================================================
calc_or_original <- function(df, region_residues, cat) {
  df[, in_region := Pos_real %in% region_residues]
  df[, is_cat := category == cat]
  
  a <- sum(df$in_region & df$is_cat)
  b <- sum(!df$in_region & df$is_cat)
  c <- sum(df$in_region & !df$is_cat)
  d <- sum(!df$in_region & !df$is_cat)
  
  tab <- matrix(c(a, b, c, d), nrow = 2)
  
  ft <- fisher.test(tab)
  
  list(OR = unname(ft$estimate), 
       p = ft$p.value, 
       OR_low = ft$conf.int[1], 
       OR_high = ft$conf.int[2])
}

# =========================================================
run_pair_all_regions <- function(pair_name, input_files, anno, structure_regions) {
  
  assays <- strsplit(pair_name, " vs ")[[1]]
  x <- assays[1]
  y <- assays[2]
  
  prepared <- prepare_merged_data_with_FDR(
    input_files[[x]], input_files[[y]], x, y, anno
  )
  
  df <- prepared$data
  threshold_x <- prepared$threshold_x
  threshold_y <- prepared$threshold_y
  
  # 
  df <- classify_two_step(df, threshold_x, threshold_y, x, y)
  
  # 
  df[, category := fifelse(
    final_classification %in% c("Both promoting", "Both disrupting"),
    "Correlated",
    fifelse(
      final_classification %in% c("Promoting in X / Disrupting in Y",
                                  "Disrupting in X / Promoting in Y"),
      "Anti-correlated",
      "Other"
    )
  )]
  
  # 
  all_results <- list()
  
  for(region_name in names(structure_regions)) {
    region_residues <- structure_regions[[region_name]]
    
    # 
    plot_df <- df[, .(
      frac = mean(Pos_real %in% region_residues),
      n = .N,
      se = sqrt(mean(Pos_real %in% region_residues) * (1 - mean(Pos_real %in% region_residues)) / .N)
    ), by = category]
    plot_df[, region := region_name]
    plot_df[, pair := pair_name]
    
    # 
    categories_to_test <- c("Correlated", "Anti-correlated", "Other")
    or_list <- lapply(categories_to_test, function(cat) {
      res <- calc_or_original(df, region_residues, cat)
      data.table(
        pair = pair_name,
        region = region_name,
        category = cat,
        OR = res$OR,
        OR_low = res$OR_low,
        OR_high = res$OR_high,
        p = res$p
      )
    })
    or_df <- rbindlist(or_list)
    
    all_results[[region_name]] <- list(plot = plot_df, or = or_df)
  }
  
  # 
  combined_plot <- rbindlist(lapply(all_results, `[[`, "plot"))
  combined_or <- rbindlist(lapply(all_results, `[[`, "or"))
  
  # 
  cat("\n", pair_name, "\n")
  cat("  Total mutations:", nrow(df), "\n")
  cat("  Correlated:", sum(df$category == "Correlated"), "\n")
  cat("  Anti-correlated:", sum(df$category == "Anti-correlated"), "\n")
  cat("  Other:", sum(df$category == "Other"), "\n")
  
  list(plot = combined_plot, or = combined_or, full_data = df)
}

# =========================================================
input_files <- list(
  RAF1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
  K55  = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K55.txt",
  K27  = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K27.txt",
  K13  = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  K19  = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt"
)

anno <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"

# 
structure_regions <- list(
  "Core" = core_residues,
  "Surface" = surface_residues,
  "NBP" = NBP,
  "Functional Loop" = functional_loop_residues,
  "Beta Sheets" = beta_sheets,
  "Alpha Helices" = alpha_helices
)

# =========================================================
target_pair <- "K27 vs K13"

res <- run_pair_all_regions(target_pair, input_files, anno, structure_regions)

plot_df <- res$plot
or_df <- res$or

# =========================================================

# 
categories_to_plot <- c("Correlated", "Anti-correlated", "Other")

#
color_map <- c(
  "Correlated" = "#F4AD0C",
  "Anti-correlated" = "#1B38A6",
  "Other" = "grey80"
)

#
plot_subset <- plot_df[category %in% categories_to_plot]
plot_subset[, region := factor(region, levels = names(structure_regions))]
plot_subset[, category := factor(category, levels = categories_to_plot)]

#
or_subset <- or_df[pair == target_pair & category %in% categories_to_plot]
or_subset[, region := factor(region, levels = names(structure_regions))]

#
or_subset[, x_offset := dplyr::case_when(
  category == "Correlated" ~ -0.25,
  category == "Anti-correlated" ~ 0,
  category == "Other" ~ 0.25
)]

#
fixed_y_position <- 1.2

#
max_y <- max(plot_subset$frac + plot_subset$se, na.rm = TRUE)
if (fixed_y_position > max_y) {
  max_y <- fixed_y_position + 0.1
} else {
  max_y <- max_y * 1.15
}
max_y <- max(max_y, 1.0)  

# 
or_subset[, label := paste0("OR = ", round(OR, 2), 
                            ifelse(p < 0.05, "*", ""))]

# 
p <- ggplot(plot_subset, aes(x = region, y = frac, fill = category)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  
  # 
  geom_errorbar(aes(ymin = frac - se, ymax = frac + se),
                position = position_dodge(0.8), 
                width = 0.15, 
                size = 0.8,
                color = "#F1DD10",
                alpha = 0.8,
                na.rm = TRUE) +
  
  # 
  geom_errorbar(aes(ymin = frac - se, ymax = frac - se),
                position = position_dodge(0.8),
                width = 0.3,
                size = 0.8,
                color = "#F1DD10",
                alpha = 0.8,
                na.rm = TRUE) +
  geom_errorbar(aes(ymin = frac + se, ymax = frac + se),
                position = position_dodge(0.8),
                width = 0.3,
                size = 0.8,
                color = "#F1DD10",
                alpha = 0.8,
                na.rm = TRUE) +
  
  # 
  geom_text(data = or_subset, 
            aes(x = as.numeric(region) + x_offset, 
                y = fixed_y_position, 
                label = label),
            size = 3.5,
            angle = 45,
            hjust = 0.5,
            vjust = 0,
            color = ifelse(or_subset$category == "Correlated", "#F4AD0C",
                           ifelse(or_subset$category == "Anti-correlated", "#1B38A6", "grey50"))) +
  
  scale_fill_manual(values = color_map) +
  
  scale_y_continuous(limits = c(0, max_y), 
                     expand = expansion(mult = c(0, 0.02)),
                     breaks = seq(0, 1, 0.2)) +
  
  # 
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.5, size = 0.8) +
  
  labs(y = "Fraction of mutations in region",
       x = "Structural region",
       title = "K27 vs K13") +
  
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 15),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 40, r = 10, b = 10, l = 10),
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# 
print(p)

# =========================================================
output_dir <- "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result5_figure6/"

ggsave(paste0(output_dir, "barplot_K27_vs_K13_core_surface_nbp_fl_bsheet_ahelix_enrichment.pdf"),
       plot = p,
       width = 10,
       height = 6,
       device = cairo_pdf)

# =========================================================
or_output <- or_df[, .(pair, region, category, OR, OR_low, OR_high, p)]
or_output[, p_signif := ifelse(p < 0.05, "*", "ns")]
or_output[, OR_text := paste0(round(OR, 2), ifelse(p < 0.05, "*", ""))]

fwrite(or_output, 
       paste0(output_dir, "OR_values_K27_vs_K13_all_six_regions.csv"))

# 
frac_output <- plot_df[, .(pair, region, category, frac, se, n)]
fwrite(frac_output,
       paste0(output_dir, "fraction_summary_K27_vs_K13_all_six_regions.csv"))

cat("\n\nAnalysis completed for", target_pair, "\n")
cat("Results saved to:", output_dir, "\n")
cat("\nAnalyzed regions:\n")
for(region in names(structure_regions)) {
  cat("  -", region, ":", length(structure_regions[[region]]), "residues\n")
}
