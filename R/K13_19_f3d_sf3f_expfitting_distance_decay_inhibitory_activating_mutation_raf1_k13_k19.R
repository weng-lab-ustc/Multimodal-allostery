library(data.table)
library(ggplot2)
library(dplyr)
library(grid)
library(krasddpcams)

# ============================================
# Main plotting function 
# ============================================
plot_energy_distance_decay_directional_noFilter <- function(
    input,
    assay_sele,
    anno_file,
    x_range = c(0, 35),
    y_range = c(-1.5, 3)
){
  
  # ============================================
  # 1. Read data
  # ============================================
  data <- fread(input)
  data[, Pos_real := Pos + 1]
  
  anno <- fread(anno_file)
  anno[, Pos_real := Pos]
  
  anno_final <- merge(anno, data, by = "Pos_real", all = FALSE)
  
  # Distance column
  x_col <- paste0("scHAmin_ligand_", assay_sele)
  
  # ============================================
  # 2. Define directional groups by ΔΔG sign 
  # ============================================
  
  # inhibit binding: ΔΔG > 0 (use ALL mutations, use absolute value)
  df_inhibit <- anno_final[
    `mean_kcal/mol` > 0,
    .(
      x = get(x_col),
      y = abs(`mean_kcal/mol`)
    )
  ]
  
  # activate binding: ΔΔG < 0 (use ALL mutations, use original negative values)
  df_activate <- anno_final[
    `mean_kcal/mol` < 0,
    .(
      x = get(x_col),
      y = `mean_kcal/mol`  # Keep original negative values
    )
  ]
  
  # remove NA
  df_inhibit <- df_inhibit[complete.cases(df_inhibit)]
  df_activate <- df_activate[complete.cases(df_activate)]
  
  # remove interface residues (<5 Å) for fitting
  df_inhibit_fit <- df_inhibit[x >= 5]
  df_activate_fit <- df_activate[x >= 5]
  
  # ============================================
  # 3. Exponential fitting: inhibit (positive values, decay to 0)
  # ============================================
  fit_inhibit <- tryCatch(
    nls(
      y ~ a * exp(b * x),
      data = df_inhibit_fit,
      start = list(a = 1, b = -0.1)
    ),
    error = function(e) NULL
  )
  
  fit_inhibit_df <- data.frame()
  annotation_inhibit <- NULL
  
  if(!is.null(fit_inhibit)){
    
    x_seq <- seq(
      min(df_inhibit_fit$x),
      max(df_inhibit_fit$x),
      length.out = 200
    )
    
    fit_inhibit_df <- data.frame(
      x = x_seq,
      y = predict(
        fit_inhibit,
        newdata = data.frame(x = x_seq)
      )
    )
    
    coefs <- summary(fit_inhibit)$coefficients
    
    annotation_inhibit <- paste0(
      "Inhibit binding\n",
      "a = ", round(coefs["a", "Estimate"], 3),
      "\n",
      "b = ", round(coefs["b", "Estimate"], 3)
    )
  }
  
  # ============================================
  # 4. Exponential fitting: activate (negative values, decay to 0 from below)
  #    Formula: y ~ a * exp(b * x), where b should be positive (since negative values approach 0)
  # ============================================
  fit_activate <- tryCatch(
    nls(
      y ~ a * exp(b * x),
      data = df_activate_fit,
      start = list(a = -1, b = -0.1)  # Start with negative a
    ),
    error = function(e) NULL
  )
  
  fit_activate_df <- data.frame()
  annotation_activate <- NULL
  
  if(!is.null(fit_activate)){
    
    x_seq <- seq(
      min(df_activate_fit$x),
      max(df_activate_fit$x),
      length.out = 200
    )
    
    fit_activate_df <- data.frame(
      x = x_seq,
      y = predict(
        fit_activate,
        newdata = data.frame(x = x_seq)
      )
    )
    
    coefs <- summary(fit_activate)$coefficients
    
    annotation_activate <- paste0(
      "Stabilize binding\n",
      "a = ", round(coefs["a", "Estimate"], 3),
      "\n",
      "b = ", round(coefs["b", "Estimate"], 3)
    )
  }
  
  # ============================================
  # 5. Plot
  # ============================================
  p <- ggplot() +
    
    # inhibit mutations (ALL, positive values)
    geom_point(
      data = df_inhibit,
      aes(x = x, y = y),
      color = "#F4270C",
      alpha = 0.25,
      size = 1.5
    ) +
    
    # activate mutations (ALL, negative values)
    geom_point(
      data = df_activate,
      aes(x = x, y = y),
      color = "#1B38A6",
      alpha = 0.25,
      size = 1.5
    ) +
    
    # inhibit fit (positive curve)
    geom_line(
      data = fit_inhibit_df,
      aes(x = x, y = y),
      color = "gray40",
      linewidth = 1
    ) +
    
    # activate fit (negative curve)
    geom_line(
      data = fit_activate_df,
      aes(x = x, y = y),
      color = "gray40",
      linewidth = 1
    ) +
    
    geom_vline(
      xintercept = 5,
      linetype = "dashed",
      color = "gray50"
    ) +
    
    scale_x_continuous(
      limits = x_range,
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      limits = y_range,
      expand = c(0, 0)
    ) +
    
    theme_classic(base_size = 10) +
    
    labs(
      x = paste0("Distance to ", assay_sele, " (Å)"),
      y = paste0("Binding ΔΔG (", assay_sele, ") (kcal/mol)")
    ) +
    
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10)
    )
  
  # ============================================
  # 6. Add annotations
  # ============================================
  if(!is.null(annotation_inhibit)){
    
    p <- p + annotate(
      "text",
      x = max(x_range) * 0.95,
      y = max(y_range) * 0.95,
      label = annotation_inhibit,
      hjust = 1,
      vjust = 1,
      size = 2.8,
      color = "#F4270C"
    )
  }
  
  if(!is.null(annotation_activate)){
    
    p <- p + annotate(
      "text",
      x = max(x_range) * 0.95,
      y = max(y_range) * 0.65,
      label = annotation_activate,
      hjust = 1,
      vjust = 1,
      size = 2.8,
      color = "#1B38A6"
    )
  }
  
  return(p)
}

# ============================================
# Run RAF1 
# ============================================
plot_raf1 <- plot_energy_distance_decay_directional_noFilter(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
  assay_sele = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
)

print(plot_raf1)

ggsave(
  "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3_4/RAF1_directional_decay_noFDR v3.pdf",
  plot_raf1,
  device = cairo_pdf,
  width = 2.5,
  height = 2.5,
  units = "in",
  dpi = 300
)

# ============================================
# Run K13 
# ============================================
plot_K13 <- plot_energy_distance_decay_directional_noFilter(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  assay_sele = "K13",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
)

print(plot_K13)

ggsave(
  "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3_4/K13_directional_decay_noFDR v3.pdf",
  plot_K13,
  device = cairo_pdf,
  width = 2.5,
  height = 2.5,
  units = "in",
  dpi = 300
)

# ============================================
# Run K19 
# ============================================
plot_K19 <- plot_energy_distance_decay_directional_noFilter(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt",
  assay_sele = "K19",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
)

print(plot_K19)

ggsave(
  "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3_4/K19_directional_decay_noFDR v3.pdf",
  plot_K19,
  device = cairo_pdf,
  width = 2.5,
  height = 2.5,
  units = "in",
  dpi = 300
)
