library(data.table)
library(ggplot2)
library(dplyr)
library(grid)

plot_energy_distance_decay_expfit_directional_noFDR <- function(
    input,
    assay_sele,
    contact_shell,
    x_range = c(0, 10),
    y_range = c(-1.5, 3)  
) {
  
  # ===============================
  # 1. Read data
  # ===============================
  data <- fread(input)
  data <- data[, Pos_real := Pos + 1]
  data <- data[, c(20:23)]
  colnames(data)[1:3] <- paste0(colnames(data)[1:3], "_", assay_sele)
  
  # ===============================
  # 2. Merge with contact shell
  # ===============================
  anno_final <- merge(contact_shell, data, by = "Pos_real", all = FALSE)
  
  # ===============================
  # 3. Select contact shell distance
  # ===============================
  x_col <- paste0(assay_sele, "_contact_shell")
  mean_col <- paste0("mean_kcal/mol_", assay_sele)
  
  # ===============================
  # 4. Define directional groups (ALL mutations, NO p-value filtering)
  # ===============================
  # Destabilize binding: ΔΔG > 0 (use ALL mutations, use absolute value)
  df_inhibit <- anno_final[
    get(mean_col) > 0,
    .(
      x = get(x_col),
      y = abs(get(mean_col))  
    )
  ]
  
  # Stabilize binding: ΔΔG < 0 (use ALL mutations, use original negative values)
  df_stabilize <- anno_final[
    get(mean_col) < 0,
    .(
      x = get(x_col),
      y = get(mean_col)  
    )
  ]
  
  # Remove NA/infinite values and x == 0 (interface)
  df_inhibit <- df_inhibit[complete.cases(df_inhibit) & is.finite(x) & is.finite(y) & x > 0]
  df_stabilize <- df_stabilize[complete.cases(df_stabilize) & is.finite(x) & is.finite(y) & x > 0]
  
  # Separate interface (x == 1) from shell (x > 1) for fitting
  df_inhibit_fit <- df_inhibit[x > 1]
  df_stabilize_fit <- df_stabilize[x > 1]
  
  # ===============================
  # 5. Exponential fitting: destabilize (ΔΔG > 0, positive values)
  #    Formula: y ~ a * exp(b * x), b < 0 for decay to 0
  # ===============================
  fit_inhibit <- tryCatch(
    nls(y ~ a * exp(b * x),
        data = df_inhibit_fit,
        start = list(a = 1, b = -0.1)),
    error = function(e) NULL
  )
  
  fit_inhibit_df <- data.frame()
  annotation_inhibit <- NULL
  
  if (!is.null(fit_inhibit)) {
    x_seq <- seq(min(df_inhibit_fit$x), max(df_inhibit_fit$x), length.out = 200)
    fit_inhibit_df <- data.frame(
      x = x_seq,
      y = predict(fit_inhibit, newdata = data.frame(x = x_seq))
    )
    
    coefs <- summary(fit_inhibit)$coefficients
    
    annotation_inhibit <- paste0(
      "Destabilize (ΔΔG>0):\n",
      "a = ", round(coefs["a", "Estimate"], 3),
      "\nb = ", round(coefs["b", "Estimate"], 3)
    )
  }
  
  # ===============================
  # 6. Exponential fitting: stabilize (ΔΔG < 0, negative values)
  #    Formula: y ~ a * exp(b * x), where a < 0, b < 0 for decay from negative to 0
  # ===============================
  fit_stabilize <- tryCatch(
    nls(y ~ a * exp(b * x),
        data = df_stabilize_fit,
        start = list(a = -1, b = -0.1)),  # Start with negative a
    error = function(e) NULL
  )
  
  fit_stabilize_df <- data.frame()
  annotation_stabilize <- NULL
  
  if (!is.null(fit_stabilize)) {
    x_seq <- seq(min(df_stabilize_fit$x), max(df_stabilize_fit$x), length.out = 200)
    fit_stabilize_df <- data.frame(
      x = x_seq,
      y = predict(fit_stabilize, newdata = data.frame(x = x_seq))
    )
    
    coefs <- summary(fit_stabilize)$coefficients
    
    annotation_stabilize <- paste0(
      "Stabilize (ΔΔG<0):\n",
      "a = ", round(coefs["a", "Estimate"], 3),
      "\nb = ", round(coefs["b", "Estimate"], 3)
    )
  }
  
  # ===============================
  # 7. Calculate median per distance for points
  # ===============================
  # For destabilize mutations - median per distance (positive values)
  df_inhibit_median <- df_inhibit %>%
    group_by(x) %>%
    summarise(y_median = median(y, na.rm = TRUE), .groups = "drop")
  
  # For stabilize mutations - median per distance (negative values)
  df_stabilize_median <- df_stabilize %>%
    group_by(x) %>%
    summarise(y_median = median(y, na.rm = TRUE), .groups = "drop")
  
  # ===============================
  # 8. Plot
  # ===============================
  p <- ggplot() +
    
    # Individual points - destabilize (ΔΔG > 0, positive)
    geom_point(data = df_inhibit, aes(x = x, y = y),
               alpha = 0.15, size = 1.5, color = "#F4270C") +
    
    # Individual points - stabilize (ΔΔG < 0, negative)
    geom_point(data = df_stabilize, aes(x = x, y = y),
               alpha = 0.15, size = 1.5, color = "#1B38A6") +
    
    # Median points - destabilize
    geom_point(data = df_inhibit_median, aes(x = x, y = y_median),
               color = "#F4270C", size = 2, shape = 16) +
    
    # Median points - stabilize
    geom_point(data = df_stabilize_median, aes(x = x, y = y_median),
               color = "#1B38A6", size = 2, shape = 16) +
    
    # Interface boundary
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    
    # Horizontal line at y=0
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.3) +
    
    # Exponential fit - destabilize (positive curve)
    geom_line(data = fit_inhibit_df, aes(x = x, y = y),
              color = "gray40", linewidth = 1) +
    
    # Exponential fit - stabilize (negative curve)
    geom_line(data = fit_stabilize_df, aes(x = x, y = y),
              color = "gray40", linewidth = 1) +
    
    # Axis settings
    scale_x_continuous(limits = x_range, expand = c(0, 0), 
                       breaks = seq(0, 10, by = 2)) +
    scale_y_continuous(limits = y_range, expand = c(0, 0), 
                       breaks = seq(-3, 3, by = 1)) +
    
    theme_classic(base_size = 10) +
    
    labs(
      x = paste0("Contact shell distance to ", assay_sele, " (Å)"),
      y = paste0("Binding ΔΔG (", assay_sele, ") (kcal/mol)")  
    ) +
    
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10)
    )
  
  # ===============================
  # 9. Add annotations
  # ===============================
  if (!is.null(annotation_inhibit)) {
    p <- p + annotate("text",
                      x = max(x_range) * 0.95,
                      y = max(y_range) * 0.95,
                      label = annotation_inhibit,
                      hjust = 1,
                      vjust = 1,
                      size = 2.5,
                      color = "#F4270C")
  }
  
  if (!is.null(annotation_stabilize)) {
    p <- p + annotate("text",
                      x = max(x_range) * 0.95,
                      y = max(y_range) * 0.65,
                      label = annotation_stabilize,
                      hjust = 1,
                      vjust = 1,
                      size = 2.5,
                      color = "#1B38A6")
  }
  
  return(p)
}

# ===============================
# Load contact shell data
# ===============================
contact_shell <- fread("C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/distance_contacts_analysis/5binder_contact_shell2.csv")

# ===============================
# Run RAF1
# ===============================
plot_raf1 <- plot_energy_distance_decay_expfit_directional_noFDR(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
  assay_sele = "RAF1",
  contact_shell = contact_shell,
  y_range = c(-1.5, 3)  
)

print(plot_raf1)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3_4/RAF1_contact_shell_decay_directional_noFDR v2.pdf", 
       plot_raf1, device = cairo_pdf, height = 2.5, width = 2.5)

# ===============================
# Run K13
# ===============================
plot_k13 <- plot_energy_distance_decay_expfit_directional_noFDR(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  assay_sele = "K13",
  contact_shell = contact_shell,
  y_range = c(-1.5, 3)
)

print(plot_k13)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3_4/K13_contact_shell_decay_directional_noFDR v2.pdf", 
       plot_k13, device = cairo_pdf, height = 2.5, width = 2.5)

# ===============================
# Run K19
# ===============================
plot_k19 <- plot_energy_distance_decay_expfit_directional_noFDR(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt",
  assay_sele = "K19",
  contact_shell = contact_shell,
  y_range = c(-1.5, 3)
)

print(plot_k19)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3_4/K19_contact_shell_decay_directional_noFDR v2.pdf", 
       plot_k19, device = cairo_pdf, height = 2.5, width = 2.5)
