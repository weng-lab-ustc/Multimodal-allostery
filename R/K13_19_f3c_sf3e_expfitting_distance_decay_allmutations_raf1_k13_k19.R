library(data.table)
library(ggplot2)
library(dplyr)
library(grid)


plot_energy_distance_decay_expfit_self <- function(input, assay_sele, anno_file,
                                                   x_range = c(0, 35), y_range = c(0, 3),
                                                   plot_width = 4, plot_height = 4) {
  
  # ===============================
  # 1. Read data
  # ===============================
  data <- fread(input)
  data <- data[, Pos_real := Pos + 1]
  data <- data[, c(20:23)]
  colnames(data)[1:3] <- paste0(colnames(data)[1:3], "_", assay_sele)
  
  anno <- fread(anno_file)
  anno <- anno[, Pos_real := Pos]
  
  anno_final <- merge(anno, data, by = "Pos_real", all = FALSE)
  
  # ===============================
  # 2. Select distance + energy
  # ===============================
  x_col <- paste0("scHAmin_ligand_", assay_sele)
  y_col <- paste0("mean_kcal/mol_", assay_sele)
  title <- paste0("Distance to ", assay_sele, " (Å)")
  
  xvector <- anno_final[[x_col]]
  yvector <- abs(anno_final[[y_col]])
  
  df <- data.frame(x = xvector, y = yvector)
  df <- df[complete.cases(df), ]
  df <- df[df$x > 0, ]
  
  # ===============================
  # 3. Interface split
  # ===============================
  df_bi <- df[df$x < 5, ]
  df_non_bi <- df[df$x >= 5, ]
  
  # ===============================
  # 4. Residue-level aggregation
  # ===============================
  df_residue <- df %>%
    group_by(x) %>%
    summarise(y = median(y, na.rm = TRUE), .groups = "drop")
  
  df_residue_non_bi <- df_residue[df_residue$x >= 5, ]
  
  # ===============================
  # 5. Mutation-level fitting
  # ===============================
  fit_mut <- tryCatch(
    nls(y ~ a * exp(b * x),
        data = df_non_bi,
        start = list(a = 1, b = -0.1)),
    error = function(e) NULL
  )
  
  fit_mut_df <- data.frame()
  annotation_mut <- NULL
  
  if (!is.null(fit_mut)) {
    
    x_seq <- seq(min(df_non_bi$x), max(df_non_bi$x), length.out = 200)
    
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
  # 6. Residue-level fitting
  # ===============================
  fit_res <- tryCatch(
    nls(y ~ a * exp(b * x),
        data = df_residue_non_bi,
        start = list(a = 1, b = -0.1)),
    error = function(e) NULL
  )
  
  fit_res_df <- data.frame()
  annotation_res <- NULL
  
  if (!is.null(fit_res)) {
    
    x_seq <- seq(min(df_residue_non_bi$x), max(df_residue_non_bi$x), length.out = 200)
    
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
  # 7. Median points
  # ===============================
  df_median_non_bi <- df_non_bi %>%
    group_by(x) %>%
    summarise(y = median(y, na.rm = TRUE), .groups = "drop")
  
  df_median_bi <- df_bi %>%
    group_by(x) %>%
    summarise(y = median(y, na.rm = TRUE), .groups = "drop")
  
  # ===============================
  # 8. Plot with adjustable text size
  # ===============================
  p <- ggplot() +
    
    geom_point(data = df_non_bi, aes(x = x, y = y),
               alpha = 0.1, size = 1.5, color = "#75C2F6") +  
    
    geom_point(data = df_bi, aes(x = x, y = y),
               alpha = 0.1, size = 1.5, color = "#FFB6C1") +  
    
    geom_point(data = df_median_non_bi, aes(x = x, y = y),
               color = "#1B38A6", size = 2) +
    
    geom_point(data = df_median_bi, aes(x = x, y = y),
               color = "#8B0000", size = 2) +
    
    geom_line(data = fit_mut_df,
              aes(x = x, y = y),
              color = "#75C2F6", linewidth = 0.8) +  # 稍细的线
    
    geom_line(data = fit_res_df,
              aes(x = x, y = y),
              color = "#1B38A6", linewidth = 0.8, linetype = "dashed") +
    
    geom_vline(xintercept = 5, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    
    scale_x_continuous(limits = x_range, expand = c(0, 0)) +
    scale_y_continuous(limits = y_range, expand = c(0, 0)) +
    
    theme_classic(base_size = 10) + 
    
    labs(
      x = title,
      y = paste0("Binding |ΔΔG| (", assay_sele, ") (kcal/mol)")
    ) +
    
    theme(
      axis.title = element_text(size = 10),   
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 10)
    )
  
  # ===============================
  # 9. Add colored annotations (adjusted position)
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
# Run RAF1
# ===============================
plot_raf1 <- plot_energy_distance_decay_expfit_self(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
  assay_sele = "RAF1",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
)

print(plot_raf1)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3/RAF1_distance_decay all mutations.pdf", 
       plot_raf1 , device = cairo_pdf,
       width = 2.5,   
       height = 2.5,  
       units = "in",
       dpi = 300)


# ===============================
# Run K13
# ===============================
plot_k13 <- plot_energy_distance_decay_expfit_self(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  assay_sele = "K13",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
)

print(plot_k13)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3/K13_distance_decay all mutations.pdf", 
       plot_k13 , device = cairo_pdf,
       width = 2.5,   
       height = 2.5,  
       units = "in",
       dpi = 300)


# ===============================
# Run K19
# ===============================
plot_k19 <- plot_energy_distance_decay_expfit_self(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt",
  assay_sele = "K19",
  anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"
)

print(plot_k19)


ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3/K19_distance_decay all mutations.pdf", 
       plot_k19 , device = cairo_pdf,
       width = 2.5,  
       height = 2.5,  
       units = "in",
       dpi = 300)

