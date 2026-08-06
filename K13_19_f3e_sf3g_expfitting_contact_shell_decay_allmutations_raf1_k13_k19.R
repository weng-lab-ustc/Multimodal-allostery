library(data.table)
library(ggplot2)
library(dplyr)
library(grid)

plot_energy_distance_decay_expfit_self <- function(input, assay_sele, contact_shell,
                                                   x_range = c(0, 10), y_range = c(0, 3)) {
  
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
  # 3. Select contact shell distance + energy
  # ===============================
  x_col <- paste0(assay_sele, "_contact_shell")
  y_col <- paste0("mean_kcal/mol_", assay_sele)
  title <- paste0("Contact shell to ", assay_sele)
  
  xvector <- anno_final[[x_col]]
  yvector <- abs(anno_final[[y_col]])
  
  df <- data.frame(Pos_real = anno_final$Pos_real, x = xvector, y = yvector)
  df <- df[complete.cases(df), ]
  df <- df[df$x > 0, ]
  
  # ===============================
  # 4. Calculate median per residue (per position)
  # ===============================
  # For each residue position, calculate the median energy across all mutations
  df_per_residue <- df %>%
    group_by(Pos_real, x) %>%
    summarise(y_median = median(y, na.rm = TRUE),
              .groups = "drop")
  
  # ===============================
  # 5. Split interface (x == 1) and non-interface (x > 1)
  # ===============================
  df_interface <- df_per_residue[df_per_residue$x == 1, ]  # x = 1 is interface
  df_shell <- df_per_residue[df_per_residue$x > 1, ]       # x > 1 is shell
  
  # ===============================
  # 6. Mutation-level fitting (using all individual points, x > 1 only)
  # ===============================
  df_mut_fit <- df[df$x > 1, ]
  
  fit_mut <- tryCatch(
    nls(y ~ a * exp(b * x),
        data = df_mut_fit,
        start = list(a = 1, b = -0.1)),
    error = function(e) NULL
  )
  
  fit_mut_df <- data.frame()
  annotation_mut <- NULL
  
  if (!is.null(fit_mut)) {
    
    x_seq <- seq(min(df_mut_fit$x), max(df_mut_fit$x), length.out = 200)
    
    fit_mut_df <- data.frame(
      x = x_seq,
      y = predict(fit_mut, newdata = data.frame(x = x_seq))
    )
    
    coefs <- summary(fit_mut)$coefficients
    p_val_b <- coefs["b", "Pr(>|t|)"]
    
    annotation_mut <- paste0(
      "Mutation-level:\n",
      "a = ", round(coefs["a","Estimate"], 3),
      "\nb = ", round(coefs["b","Estimate"], 3),
      "\np = ", signif(p_val_b, 3)
    )
  }
  
  # ===============================
  # 7. Residue-level fitting (using per-residue medians, x > 1 only)
  # ===============================
  df_res_fit <- df_shell[, c("x", "y_median")]
  colnames(df_res_fit)[2] <- "y"
  
  fit_res <- tryCatch(
    nls(y ~ a * exp(b * x),
        data = df_res_fit,
        start = list(a = 1, b = -0.1)),
    error = function(e) NULL
  )
  
  fit_res_df <- data.frame()
  annotation_res <- NULL
  
  if (!is.null(fit_res)) {
    
    x_seq <- seq(min(df_res_fit$x), max(df_res_fit$x), length.out = 200)
    
    fit_res_df <- data.frame(
      x = x_seq,
      y = predict(fit_res, newdata = data.frame(x = x_seq))
    )
    
    coefs <- summary(fit_res)$coefficients
    p_val_b <- coefs["b", "Pr(>|t|)"]
    
    annotation_res <- paste0(
      "Residue-level:\n",
      "a = ", round(coefs["a","Estimate"], 3),
      "\nb = ", round(coefs["b","Estimate"], 3),
      "\np = ", signif(p_val_b, 3)
    )
  }
  
  # ===============================
  # 8. Prepare median points for plotting
  # ===============================
  df_median_shell <- df_shell
  colnames(df_median_shell)[3] <- "y"
  
  df_median_interface <- df_interface
  colnames(df_median_interface)[3] <- "y"
  
  # ===============================
  # 9. Plot
  # ===============================
  p <- ggplot() +
    
    # Individual mutation points (shell region, x > 1)
    geom_point(data = df[df$x > 1, ], aes(x = x, y = y),
               alpha = 0.1, size = 2, color = "#75C2F6") +
    
    # Individual mutation points (interface, x == 1)
    geom_point(data = df[df$x == 1, ], aes(x = x, y = y),
               alpha = 0.1, size = 2, color = "#FFB6C1") +
    
    # Median points for shell (x > 1)
    geom_point(data = df_median_shell, aes(x = x, y = y),
               color = "#1B38A6", size = 2) +
    
    # Median points for interface (x == 1)
    geom_point(data = df_median_interface, aes(x = x, y = y),
               color = "#8B0000", size = 2) +
    
    # Mutation-level fit line (using all individual points)
    geom_line(data = fit_mut_df,
              aes(x = x, y = y),
              color = "#75C2F6", linewidth = 1) +
    
    # Residue-level fit line (using medians)
    geom_line(data = fit_res_df,
              aes(x = x, y = y),
              color = "#1B38A6", linewidth = 1, linetype = "dashed") +
    
    # Axis settings - with integer breaks
    scale_x_continuous(limits = x_range, expand = c(0, 0), breaks = seq(0, 10, by = 1)) +
    scale_y_continuous(limits = y_range, expand = c(0, 0), breaks = seq(0, 3, by = 1)) +
    
    theme_classic(base_size = 10) +
    
    labs(
      x = title,
      y = paste0("Binding |ΔΔG| (", assay_sele, ") (kcal/mol)")
    )
  
  # ===============================
  # 10. Add annotations
  # ===============================
  if (!is.null(annotation_mut)) {
    p <- p + annotate("text",
                      x = max(x_range)*0.95,
                      y = max(y_range)*0.95,
                      label = annotation_mut,
                      hjust = 1,
                      vjust = 1,
                      size = 2.8,
                      color = "#75C2F6")
  }
  
  if (!is.null(annotation_res)) {
    p <- p + annotate("text",
                      x = max(x_range)*0.95,
                      y = max(y_range)*0.60,
                      label = annotation_res,
                      hjust = 1,
                      vjust = 1,
                      size = 2.8,
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
plot_raf1 <- plot_energy_distance_decay_expfit_self(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
  assay_sele = "RAF1",
  contact_shell = contact_shell
)

print(plot_raf1)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3/RAF1_contact_shell_decay all mutations.pdf", 
       plot_raf1, device = cairo_pdf, height = 2.5, width = 2.5)

# ===============================
# Run K13
# ===============================
plot_k13 <- plot_energy_distance_decay_expfit_self(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  assay_sele = "K13",
  contact_shell = contact_shell
)

print(plot_k13)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3/K13_contact_shell_decay all mutations.pdf", 
       plot_k13, device = cairo_pdf, height = 2.5, width = 2.5)

# ===============================
# Run K19
# ===============================
plot_k19 <- plot_energy_distance_decay_expfit_self(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt",
  assay_sele = "K19",
  contact_shell = contact_shell
)

print(plot_k19)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3/K19_contact_shell_decay all mutations.pdf", 
       plot_k19, device = cairo_pdf, height = 2.5, width = 2.5)
