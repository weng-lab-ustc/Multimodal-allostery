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
# STEP 7 — RAF1 vs K27 anticorrelated mutations
# =========================================================

get_RAF1_vs_K27_anticorrelated <- function() {
  
  input_files <- list(
    RAF1 =
      "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
    
    K27 =
      "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K27.txt"
  )
  
  anno <-
    "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"
  
  prepared <- prepare_merged_data_with_FDR(
    input_files$RAF1,
    input_files$K27,
    "RAF1",
    "K27",
    anno
  )
  
  merged_data <- prepared$data
  
  threshold_RAF1 <- prepared$threshold_x
  threshold_K27  <- prepared$threshold_y
  
  merged_data[, pass_FDR_RAF1 := p_adj_x < 0.05]
  merged_data[, pass_FDR_K27  := p_adj_y < 0.05]
  
  merged_data[, direction_class :=
                classify_by_direction(
                  ddG_RAF1,
                  ddG_K27,
                  threshold_RAF1,
                  threshold_K27
                )]
  
  merged_data[, final_classification :=
                reclassify_by_FDR(
                  direction_class,
                  pass_FDR_RAF1,
                  pass_FDR_K27
                )]
  
  anticorrelated_muts <- merged_data[
    final_classification %in%
      c("Promoting in X / Disrupting in Y",
        "Disrupting in X / Promoting in Y")
  ]
  
  # Add NBP information
  anticorrelated_muts[, is_NBP := Pos_real %in% NBP_RESIDUES]
  
  cat("\nAnticorrelated mutations found:",
      nrow(anticorrelated_muts), "\n")
  cat("  - NBP residues:",
      sum(anticorrelated_muts$is_NBP), "\n")
  cat("  - Non-NBP residues:",
      sum(!anticorrelated_muts$is_NBP), "\n")
  
  return(anticorrelated_muts)
}

# =========================================================
# STEP 8 — Load K13 vs K19 data
# =========================================================

load_K13_vs_K19_data <- function() {
  
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
  
  return(list(
    data = prepared$data,
    threshold_K13 = prepared$threshold_x,
    threshold_K19 = prepared$threshold_y
  ))
}

# =========================================================
# STEP 9 — Plot (with shape distinction for NBP vs non-NBP)
# =========================================================

plot_mapped_anticorrelated <- function(
    k13_k19_data,
    raf1_k27_anticorrelated,
    xlim = c(-1.5, 3),
    ylim = c(-1.5, 3)
) {
  
  merged_data <- copy(k13_k19_data$data)
  
  threshold_K13 <- k13_k19_data$threshold_K13
  threshold_K19 <- k13_k19_data$threshold_K19
  
  # =========================================================
  # standardize mutation names again and add NBP info
  # =========================================================
  
  merged_data[, mt := toupper(trimws(mt))]
  raf1_k27_anticorrelated[, mt := toupper(trimws(mt))]
  
  # Add NBP information to merged_data for matching mutations
  merged_data[, is_NBP := Pos_real %in% NBP_RESIDUES]
  
  # Mark which mutations are in the target set
  merged_data[, is_target :=
                mt %in% raf1_k27_anticorrelated$mt]
  
  # Add NBP status from the anticorrelated data for target mutations
  merged_data[, target_NBP_status := FALSE]
  for(i in 1:nrow(raf1_k27_anticorrelated)) {
    merged_data[mt == raf1_k27_anticorrelated$mt[i], 
                target_NBP_status := raf1_k27_anticorrelated$is_NBP[i]]
  }
  
  cat("\nMatched mutations:",
      merged_data[is_target == TRUE, .N], "\n")
  
  print(
    merged_data[is_target == TRUE,
                .(mt, Pos_real, is_NBP, target_NBP_status, ddG_K13, ddG_K19)]
  )
  
  # Separate target mutations by NBP status
  target_nbp <- merged_data[is_target == TRUE & target_NBP_status == TRUE]
  target_non_nbp <- merged_data[is_target == TRUE & target_NBP_status == FALSE]
  
  # =========================================================
  # correlation test with significance
  # =========================================================
  
  cor_test <- cor.test(merged_data$ddG_K13, merged_data$ddG_K19)
  r_value <- round(cor_test$estimate, 3)
  p_value <- cor_test$p.value
  
  # Add significance stars
  sig_stars <- ifelse(p_value < 0.001, "***", 
                      ifelse(p_value < 0.01, "**", 
                             ifelse(p_value < 0.05, "*", "")))
  
  r_label <- paste0("R = ", r_value, sig_stars)
  
  p <- ggplot() +
    
    theme_classic(base_size = 20) +
    
    geom_vline(
      xintercept = c(-threshold_K13, threshold_K13),
      linetype = "dashed",
      color = "grey50"
    ) +
    
    geom_hline(
      yintercept = c(-threshold_K19, threshold_K19),
      linetype = "dashed",
      color = "grey50"
    ) +
    
    # Other mutations (grey circles)
    geom_point(
      data = merged_data[is_target == FALSE],
      aes(ddG_K13, ddG_K19),
      color = "grey80",
      size = 2,
      alpha = 0.4,
      shape = 16
    ) +
    
    # Target mutations - non-NBP (orange circles)
    geom_point(
      data = target_non_nbp,
      aes(ddG_K13, ddG_K19),
      color = "#F1DD10",
      shape = 16,  # circle for non-NBP
      size = 3
    ) +
    
    # Target mutations - NBP (orange triangles)
    geom_point(
      data = target_nbp,
      aes(ddG_K13, ddG_K19),
      color = "#F1DD10",
      shape = 17,  # triangle for NBP
      size = 3.5
    ) +
    
    # =========================================================
  # LABELS for NBP mutations only
  # =========================================================
  #geom_text_repel(
  #  data = target_nbp,
  #  aes(ddG_K13, ddG_K19, label = mt),
  #  color = "#1B38A6",
  #  fontface = "bold",
  #  size = 5,
  #  max.overlaps = Inf,
  #  force = 10,
  #  force_pull = 0.5,
  #  box.padding = 0.7,
  #  point.padding = 0.5,
  #  segment.color = "#1B38A6",
  #  segment.size = 0.5,
  #  min.segment.length = 0
  #) +
    
    annotate(
      "text",
      x = xlim[1],
      y = ylim[2],
      hjust = 0,
      vjust = 1,
      label = r_label,
      size = 6
    ) +
    
    labs(
      x = expression("Binding"~Delta*Delta*G~"(K13) (kcal/mol)"),
      y = expression("Binding"~Delta*Delta*G~"(K19) (kcal/mol)"),
      title = "RAF1 vs K27 anticorrelated mutations mapped onto K13 vs K19"
    ) +
    
    coord_cartesian(
      xlim = xlim,
      ylim = ylim
    ) +
    
    theme(
      panel.border =
        element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.8
        ),
      
      axis.text.x =
        element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ),
      
      plot.title =
        element_text(
          hjust = 0.5
        )
    )
  
  return(p)
}

# =========================================================
# STEP 10 — Run analysis
# =========================================================

cat("\n========================================================")
cat("\nSTEP 1 — RAF1 vs K27 anticorrelated mutations")
cat("\n========================================================\n")

raf1_k27_anticorrelated <-
  get_RAF1_vs_K27_anticorrelated()

cat("\n========================================================")
cat("\nSTEP 2 — Load K13 vs K19 data")
cat("\n========================================================\n")

k13_k19_data <- load_K13_vs_K19_data()

cat("\n========================================================")
cat("\nSTEP 3 — Generate plot")
cat("\n========================================================\n")

# 使用自定义坐标轴范围
mapping_plot <- plot_mapped_anticorrelated(
  k13_k19_data,
  raf1_k27_anticorrelated,
  xlim = c(-1.5, 3),
  ylim = c(-1.5, 3)
)

print(mapping_plot)

# =========================================================
# STEP 11 — Save figure
# =========================================================

output_plot_path <-
  "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/RAF1_K27_anticorrelated_mapped_to_K13_vs_K19_fixed_axis.pdf"

ggsave(
  output_plot_path,
  plot = mapping_plot,
  width = 4.5,
  height = 5,
  device = cairo_pdf
)

cat("\nPlot saved to:\n")
cat(output_plot_path, "\n")
