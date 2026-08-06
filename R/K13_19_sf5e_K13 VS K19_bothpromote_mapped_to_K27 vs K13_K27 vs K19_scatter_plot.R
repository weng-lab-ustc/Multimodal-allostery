# =========================================================
# STEP 0 — Load packages
# =========================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(grid)

# =========================================================
# STEP 1 — Global constants
# =========================================================

binding_sites_map <- list(
  RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
  K27  = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71),
  K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99,
          101, 102, 105, 106, 107, 129, 133, 136, 137, 138),
  K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99,
          101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
)

# NBP residues
NBP_RESIDUES <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 
                  57, 60, 61, 116, 117, 119, 120, 145, 146, 147)

# =========================================================
# STEP 2 — Load mutation data
# =========================================================

load_mutation_data <- function(input, assay_sele) {
  
  if (is.character(input)) {
    ddG <- fread(input)
  } else {
    ddG <- as.data.table(input)
  }
  
  ddG[, Pos_real := Pos_ref + 1]
  
  ddG[id != "WT", wt_codon := substr(id, 1, 1)]
  ddG[id != "WT", mt_codon := substr(id, nchar(id), nchar(id))]
  
  ddG[, mt := paste0(wt_codon, Pos_real, mt_codon)]
  
  ddG <- ddG[id != "WT"]
  
  result <- ddG[, .(
    mt,
    Pos_real,
    wt_codon,
    mt_codon,
    `mean_kcal/mol`,
    `std_kcal/mol`
  )]
  
  setnames(result, "mean_kcal/mol", "ddG")
  setnames(result, "std_kcal/mol", "ddG_std")
  
  result[, assay := assay_sele]
  
  # standardize mutation names
  result[, mt := toupper(trimws(mt))]
  
  return(result)
}

# =========================================================
# STEP 3 — Calculate threshold
# =========================================================

calculate_threshold <- function(data, assay_sele, anno) {
  
  anno_data <- fread(anno)
  
  site_ddG <- data[, .(
    mean_ddG =
      sum(abs(ddG) / ddG_std^2, na.rm = TRUE) /
      sum(1 / ddG_std^2, na.rm = TRUE)
  ), by = Pos_real]
  
  site_sigma <- data[, .(
    sigma = sqrt(
      1 / sum(1 / ddG_std^2, na.rm = TRUE)
    )
  ), by = Pos_real]
  
  site_stats <- merge(site_ddG, site_sigma, by = "Pos_real")
  
  site_stats <- merge(
    site_stats,
    anno_data,
    by.x = "Pos_real",
    by.y = "Pos",
    all.x = TRUE
  )
  
  scHAmin_col <- paste0("scHAmin_ligand_", assay_sele)
  
  binding_sites <- site_stats[
    get(scHAmin_col) < 5,
    Pos_real
  ]
  
  threshold <- site_stats[
    Pos_real %in% binding_sites,
    sum(abs(mean_ddG) / sigma^2, na.rm = TRUE) /
      sum(1 / sigma^2, na.rm = TRUE)
  ]
  
  return(threshold)
}

# =========================================================
# STEP 4 — Direction classification
# =========================================================

classify_by_direction <- function(ddG_x, ddG_y,
                                  threshold_x, threshold_y) {
  
  sig_x <- abs(ddG_x) > threshold_x
  sig_y <- abs(ddG_y) > threshold_y
  
  disrupt_x <- ddG_x > threshold_x
  promote_x <- ddG_x < -threshold_x
  
  disrupt_y <- ddG_y > threshold_y
  promote_y <- ddG_y < -threshold_y
  
  result <- rep("neutral", length(ddG_x))
  
  result[sig_x & sig_y & promote_x & promote_y] <-
    "Both promoting"
  
  result[sig_x & sig_y & disrupt_x & disrupt_y] <-
    "Both disrupting"
  
  result[sig_x & sig_y & promote_x & disrupt_y] <-
    "Promoting in X / Disrupting in Y"
  
  result[sig_x & sig_y & disrupt_x & promote_y] <-
    "Disrupting in X / Promoting in Y"
  
  result[sig_x & !sig_y] <-
    "Allosteric only in X"
  
  result[!sig_x & sig_y] <-
    "Allosteric only in Y"
  
  result[!sig_x & !sig_y] <-
    "Not significant (FDR >= 0.05)"
  
  return(result)
}

# =========================================================
# STEP 5 — FDR reclassification
# =========================================================

reclassify_by_FDR <- function(direction_class,
                              pass_FDR_x,
                              pass_FDR_y) {
  
  result <- direction_class
  
  for (i in 1:length(direction_class)) {
    
    if (direction_class[i] %in%
        c("Both promoting",
          "Both disrupting",
          "Promoting in X / Disrupting in Y",
          "Disrupting in X / Promoting in Y")) {
      
      if (!(pass_FDR_x[i] & pass_FDR_y[i])) {
        
        if (pass_FDR_x[i] & !pass_FDR_y[i]) {
          result[i] <- "Allosteric only in X"
          
        } else if (!pass_FDR_x[i] & pass_FDR_y[i]) {
          result[i] <- "Allosteric only in Y"
          
        } else {
          result[i] <- "Not significant (FDR >= 0.05)"
        }
      }
    }
    
    else if (direction_class[i] == "Allosteric only in X") {
      
      if (!pass_FDR_x[i]) {
        result[i] <- "Not significant (FDR >= 0.05)"
      }
    }
    
    else if (direction_class[i] == "Allosteric only in Y") {
      
      if (!pass_FDR_y[i]) {
        result[i] <- "Not significant (FDR >= 0.05)"
      }
    }
  }
  
  return(result)
}

# =========================================================
# STEP 6 — Prepare merged data
# =========================================================

prepare_merged_data_with_FDR <- function(input_x,
                                         input_y,
                                         assay_x,
                                         assay_y,
                                         anno) {
  
  data_x <- load_mutation_data(input_x, assay_x)
  data_y <- load_mutation_data(input_y, assay_y)
  
  threshold_x <- calculate_threshold(data_x, assay_x, anno)
  threshold_y <- calculate_threshold(data_y, assay_y, anno)
  
  if (!is.null(binding_sites_map[[assay_x]])) {
    data_x <- data_x[
      !(Pos_real %in% binding_sites_map[[assay_x]])
    ]
  }
  
  if (!is.null(binding_sites_map[[assay_y]])) {
    data_y <- data_y[
      !(Pos_real %in% binding_sites_map[[assay_y]])
    ]
  }
  
  data_x_clean <- data_x[, .(
    mt,
    Pos_real,
    ddG,
    ddG_std
  )]
  
  setnames(data_x_clean, "ddG", paste0("ddG_", assay_x))
  setnames(data_x_clean, "ddG_std", paste0("std_", assay_x))
  
  data_y_clean <- data_y[, .(
    mt,
    Pos_real,
    ddG,
    ddG_std
  )]
  
  setnames(data_y_clean, "ddG", paste0("ddG_", assay_y))
  setnames(data_y_clean, "ddG_std", paste0("std_", assay_y))
  
  merged_data <- merge(
    data_x_clean,
    data_y_clean,
    by = c("mt", "Pos_real")
  )
  
  pvalue_threshold <- function(av, se, threshold) {
    
    zscore <- (abs(av) - threshold) / se
    
    2 * pnorm(abs(zscore), lower.tail = FALSE)
  }
  
  merged_data[, p_x :=
                pvalue_threshold(
                  get(paste0("ddG_", assay_x)),
                  get(paste0("std_", assay_x)),
                  threshold_x
                )]
  
  merged_data[, p_y :=
                pvalue_threshold(
                  get(paste0("ddG_", assay_y)),
                  get(paste0("std_", assay_y)),
                  threshold_y
                )]
  
  merged_data[, p_adj_x := p.adjust(p_x, "BH")]
  merged_data[, p_adj_y := p.adjust(p_y, "BH")]
  
  return(list(
    data = merged_data,
    threshold_x = threshold_x,
    threshold_y = threshold_y
  ))
}

# =========================================================
# STEP 7 — Get K13 vs K19 both promoting mutations
# =========================================================

get_K13_vs_K19_both_promoting <- function() {
  
  input_files <- list(
    K13 =
      "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
    
    K19 =
      "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt"
  )
  
  anno <-
    "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"
  
  prepared <- prepare_merged_data_with_FDR(
    input_files$K13,
    input_files$K19,
    "K13",
    "K19",
    anno
  )
  
  merged_data <- prepared$data
  
  threshold_K13 <- prepared$threshold_x
  threshold_K19 <- prepared$threshold_y
  
  merged_data[, pass_FDR_K13 := p_adj_x < 0.05]
  merged_data[, pass_FDR_K19 := p_adj_y < 0.05]
  
  merged_data[, direction_class :=
                classify_by_direction(
                  ddG_K13,
                  ddG_K19,
                  threshold_K13,
                  threshold_K19
                )]
  
  merged_data[, final_classification :=
                reclassify_by_FDR(
                  direction_class,
                  pass_FDR_K13,
                  pass_FDR_K19
                )]
  
  both_promoting_muts <- merged_data[
    final_classification == "Both promoting"
  ]
  
  # Add NBP information
  both_promoting_muts[, is_NBP := Pos_real %in% NBP_RESIDUES]
  
  cat("\nK13 vs K19 Both promoting mutations found:",
      nrow(both_promoting_muts), "\n")
  cat("  - NBP residues:",
      sum(both_promoting_muts$is_NBP), "\n")
  cat("  - Non-NBP residues:",
      sum(!both_promoting_muts$is_NBP), "\n")
  
  print(both_promoting_muts[, .(mt, Pos_real, is_NBP, ddG_K13, ddG_K19)])
  
  return(both_promoting_muts)
}

# =========================================================
# STEP 8 — Load K27 vs K13 data
# =========================================================

load_K27_vs_K13_data <- function() {
  
  input_files <- list(
    K27 =
      "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K27.txt",
    
    K13 =
      "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt"
  )
  
  anno <-
    "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"
  
  prepared <- prepare_merged_data_with_FDR(
    input_files$K27,
    input_files$K13,
    "K27",
    "K13",
    anno
  )
  
  return(list(
    data = prepared$data,
    threshold_K27 = prepared$threshold_x,
    threshold_K13 = prepared$threshold_y
  ))
}

# =========================================================
# STEP 9 — Load K27 vs K19 data
# =========================================================

load_K27_vs_K19_data <- function() {
  
  input_files <- list(
    K27 =
      "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K27.txt",
    
    K19 =
      "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt"
  )
  
  anno <-
    "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"
  
  prepared <- prepare_merged_data_with_FDR(
    input_files$K27,
    input_files$K19,
    "K27",
    "K19",
    anno
  )
  
  return(list(
    data = prepared$data,
    threshold_K27 = prepared$threshold_x,
    threshold_K19 = prepared$threshold_y
  ))
}

# =========================================================
# STEP 10 — Plot both promoting mutations mapped onto K27 vs K13
# =========================================================

plot_mapped_both_promoting_K27_K13 <- function(
    k27_k13_data,
    both_promoting_muts,
    xlim = c(-1.5, 3),
    ylim = c(-1.5, 3)
) {
  
  merged_data <- copy(k27_k13_data$data)
  
  threshold_K27 <- k27_k13_data$threshold_K27
  threshold_K13 <- k27_k13_data$threshold_K13
  
  # standardize mutation names
  merged_data[, mt := toupper(trimws(mt))]
  both_promoting_muts[, mt := toupper(trimws(mt))]
  
  # Mark which mutations are in the target set
  merged_data[, is_target := mt %in% both_promoting_muts$mt]
  
  # Add NBP status
  merged_data[, target_NBP_status := FALSE]
  for(i in 1:nrow(both_promoting_muts)) {
    merged_data[mt == both_promoting_muts$mt[i], 
                target_NBP_status := both_promoting_muts$is_NBP[i]]
  }
  
  cat("\nMapped to K27 vs K13 - Matched mutations:",
      merged_data[is_target == TRUE, .N], "\n")
  
  if(merged_data[is_target == TRUE, .N] > 0) {
    print(merged_data[is_target == TRUE,
                      .(mt, target_NBP_status, ddG_K27, ddG_K13)])
  }
  
  # Separate target mutations by NBP status
  target_nbp <- merged_data[is_target == TRUE & target_NBP_status == TRUE]
  target_non_nbp <- merged_data[is_target == TRUE & target_NBP_status == FALSE]
  
  # Correlation test
  cor_test <- cor.test(merged_data$ddG_K27, merged_data$ddG_K13)
  r_value <- round(cor_test$estimate, 3)
  p_value <- cor_test$p.value
  sig_stars <- ifelse(p_value < 0.001, "***", 
                      ifelse(p_value < 0.01, "**", 
                             ifelse(p_value < 0.05, "*", "")))
  r_label <- paste0("R = ", r_value, sig_stars)
  
  p <- ggplot() +
    theme_classic(base_size = 20) +
    
    geom_vline(xintercept = c(-threshold_K27, threshold_K27),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = c(-threshold_K13, threshold_K13),
               linetype = "dashed", color = "grey50") +
    
    # Other mutations
    geom_point(data = merged_data[is_target == FALSE],
               aes(ddG_K27, ddG_K13),
               color = "grey80", size = 2, alpha = 0.4, shape = 16) +
    
    # Target mutations - non-NBP (light red circles)
    geom_point(data = target_non_nbp,
               aes(ddG_K27, ddG_K13),
               color = "#FFB0A5", shape = 16, size = 3) +
    
    # Target mutations - NBP (light red triangles)
    geom_point(data = target_nbp,
               aes(ddG_K27, ddG_K13),
               color = "#FFB0A5", shape = 17, size = 3.5) +
    
    #geom_text_repel(data = target_nbp,
    #                aes(ddG_K27, ddG_K13, label = mt),
    #                color = "#1B38A6", fontface = "bold", size = 5,
    #                max.overlaps = Inf, force = 10, box.padding = 0.7,
    #                point.padding = 0.5, segment.color = "#1B38A6",
    #                segment.size = 0.5) +
    
    annotate("text", x = xlim[1], y = ylim[2], hjust = 0, vjust = 1,
             label = r_label, size = 6) +
    
    labs(x = expression("Binding"~Delta*Delta*G~"(K27) (kcal/mol)"),
         y = expression("Binding"~Delta*Delta*G~"(K13) (kcal/mol)"),
         title = "K13 vs K19 both promoting mutations mapped onto K27 vs K13") +
    
    coord_cartesian(xlim = xlim, ylim = ylim) +
    
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  return(p)
}

# =========================================================
# STEP 11 — Plot both promoting mutations mapped onto K27 vs K19
# =========================================================

plot_mapped_both_promoting_K27_K19 <- function(
    k27_k19_data,
    both_promoting_muts,
    xlim = c(-1.5, 3),
    ylim = c(-1.5, 3)
) {
  
  merged_data <- copy(k27_k19_data$data)
  
  threshold_K27 <- k27_k19_data$threshold_K27
  threshold_K19 <- k27_k19_data$threshold_K19
  
  # standardize mutation names
  merged_data[, mt := toupper(trimws(mt))]
  both_promoting_muts[, mt := toupper(trimws(mt))]
  
  # Mark which mutations are in the target set
  merged_data[, is_target := mt %in% both_promoting_muts$mt]
  
  # Add NBP status
  merged_data[, target_NBP_status := FALSE]
  for(i in 1:nrow(both_promoting_muts)) {
    merged_data[mt == both_promoting_muts$mt[i], 
                target_NBP_status := both_promoting_muts$is_NBP[i]]
  }
  
  cat("\nMapped to K27 vs K19 - Matched mutations:",
      merged_data[is_target == TRUE, .N], "\n")
  
  if(merged_data[is_target == TRUE, .N] > 0) {
    print(merged_data[is_target == TRUE,
                      .(mt, target_NBP_status, ddG_K27, ddG_K19)])
  }
  
  # Separate target mutations by NBP status
  target_nbp <- merged_data[is_target == TRUE & target_NBP_status == TRUE]
  target_non_nbp <- merged_data[is_target == TRUE & target_NBP_status == FALSE]
  
  # Correlation test
  cor_test <- cor.test(merged_data$ddG_K27, merged_data$ddG_K19)
  r_value <- round(cor_test$estimate, 3)
  p_value <- cor_test$p.value
  sig_stars <- ifelse(p_value < 0.001, "***", 
                      ifelse(p_value < 0.01, "**", 
                             ifelse(p_value < 0.05, "*", "")))
  r_label <- paste0("R = ", r_value, sig_stars)
  
  p <- ggplot() +
    theme_classic(base_size = 20) +
    
    geom_vline(xintercept = c(-threshold_K27, threshold_K27),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = c(-threshold_K19, threshold_K19),
               linetype = "dashed", color = "grey50") +
    
    # Other mutations
    geom_point(data = merged_data[is_target == FALSE],
               aes(ddG_K27, ddG_K19),
               color = "grey80", size = 2, alpha = 0.4, shape = 16) +
    
    # Target mutations - non-NBP (light red circles)
    geom_point(data = target_non_nbp,
               aes(ddG_K27, ddG_K19),
               color = "#FFB0A5", shape = 16, size = 3) +
    
    # Target mutations - NBP (light red triangles)
    geom_point(data = target_nbp,
               aes(ddG_K27, ddG_K19),
               color = "#FFB0A5", shape = 17, size = 3.5) +
    
    #geom_text_repel(data = target_nbp,
    #                aes(ddG_K27, ddG_K19, label = mt),
    #                color = "#1B38A6", fontface = "bold", size = 5,
    #                max.overlaps = Inf, force = 10, box.padding = 0.7,
    #                point.padding = 0.5, segment.color = "#1B38A6",
    #                segment.size = 0.5) +
    
    annotate("text", x = xlim[1], y = ylim[2], hjust = 0, vjust = 1,
             label = r_label, size = 6) +
    
    labs(x = expression("Binding"~Delta*Delta*G~"(K27) (kcal/mol)"),
         y = expression("Binding"~Delta*Delta*G~"(K19) (kcal/mol)"),
         title = "K13 vs K19 both promoting mutations mapped onto K27 vs K19") +
    
    coord_cartesian(xlim = xlim, ylim = ylim) +
    
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  return(p)
}

# =========================================================
# STEP 12 — Run analysis
# =========================================================

cat("\n========================================================")
cat("\nSTEP 1 — Get K13 vs K19 both promoting mutations")
cat("\n========================================================\n")

both_promoting_muts <- get_K13_vs_K19_both_promoting()

cat("\n========================================================")
cat("\nSTEP 2 — Load K27 vs K13 data")
cat("\n========================================================\n")

k27_k13_data <- load_K27_vs_K13_data()

cat("\n========================================================")
cat("\nSTEP 3 — Load K27 vs K19 data")
cat("\n========================================================\n")

k27_k19_data <- load_K27_vs_K19_data()

cat("\n========================================================")
cat("\nSTEP 4 — Generate plot: K27 vs K13 (fixed axis)")
cat("\n========================================================\n")

mapping_plot_K27_K13 <- plot_mapped_both_promoting_K27_K13(
  k27_k13_data,
  both_promoting_muts,
  xlim = c(-1.5, 3),
  ylim = c(-1.5, 3)
)

print(mapping_plot_K27_K13)

cat("\n========================================================")
cat("\nSTEP 5 — Generate plot: K27 vs K19 (fixed axis)")
cat("\n========================================================\n")

mapping_plot_K27_K19 <- plot_mapped_both_promoting_K27_K19(
  k27_k19_data,
  both_promoting_muts,
  xlim = c(-1.5, 3),
  ylim = c(-1.5, 3)
)

print(mapping_plot_K27_K19)

# =========================================================
# STEP 13 — Save figures
# =========================================================

output_dir <- "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/"

ggsave(
  file.path(output_dir, "K13_K19_both_promoting_mapped_to_K27_vs_K13_fixed_axis.pdf"),
  plot = mapping_plot_K27_K13,
  width = 4.5,
  height = 5,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "K13_K19_both_promoting_mapped_to_K27_vs_K19_fixed_axis.pdf"),
  plot = mapping_plot_K27_K19,
  width = 4.5,
  height = 5,
  device = cairo_pdf
)

cat("\nPlots saved to:\n")
cat(file.path(output_dir, "K13_K19_both_promoting_mapped_to_K27_vs_K13_fixed_axis.pdf"), "\n")
cat(file.path(output_dir, "K13_K19_both_promoting_mapped_to_K27_vs_K19_fixed_axis.pdf"), "\n")

