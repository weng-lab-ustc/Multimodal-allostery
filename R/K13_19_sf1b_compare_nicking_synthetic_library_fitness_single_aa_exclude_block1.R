# Library Fitness Comparison Analysis - Overall Only --exclude block1
# Compare fitness data between two libraries (e.g., synthetic vs nicking)

library(wlab.block)
library(data.table)
library(ggplot2)
library(dplyr)


compare_fitness_libraries_singlemut_overall_no_block1 <- function(
    lib1_block2, lib1_block3,
    lib2_block2, lib2_block3,
    wt_aa,
    output_file = NULL,
    x_lab = "Library 1 fitness",
    y_lab = "Library 2 fitness", 
    main_title = "Comparison of fitness data between two libraries",
    point_alpha = 0.3,
    plot_width = 5,
    plot_height = 4
) {
  
  # Internal function: process single library data (only block2 + block3)
  process_library_data <- function(block2, block3, wt_aa, suffix) {
    nor_fit <- nor_fitness(block2 = block2, block3 = block3)
    nor_fit_single <- nor_fitness_single_mut(input = nor_fit)
    nor_fit_single <- pos_id(nor_fit_single, wt_aa)
    
    fitness_data <- nor_fit_single[, c(1, 40, 41, 46, 48, 50, 52)]
    colnames(fitness_data) <- c("block", 
                                paste0("fitness", suffix),
                                paste0("fitness_sigma", suffix),
                                "Pos", "wtcodon", "codon", "mt")
    return(fitness_data)
  }
  
  # Process both libraries
  cat("Processing library 1 data...\n")
  fitness_data_1 <- process_library_data(lib1_block2, lib1_block3, wt_aa, "1")
  
  cat("Processing library 2 data...\n")
  fitness_data_2 <- process_library_data(lib2_block2, lib2_block3, wt_aa, "2")
  
  # Merge data
  data <- merge(fitness_data_1, fitness_data_2, by = c("block", "Pos", "wtcodon", "codon", "mt"), all = FALSE)
  setDT(data)
  
  cat(paste("Total variants after merging:", nrow(data), "\n"))
  
  # ---- Overall correlation plot only ----
  complete_cases <- complete.cases(data$fitness1, data$fitness2)
  data_complete <- data[complete_cases, ]
  
  if (nrow(data_complete) < 2) {
    warning("Insufficient complete cases for plot")
    return(ggplot() + 
             labs(title = main_title, x = x_lab, y = y_lab) +
             annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") +
             theme_minimal())
  }
  
  # Pearson correlation
  cor_test <- cor.test(data_complete$fitness1, data_complete$fitness2, 
                       method = "pearson", use = "complete.obs")
  r_value <- round(cor_test$estimate, 3)
  p_value <- round(cor_test$p.value, 4)
  
  # Create overall plot
  p <- ggplot(data_complete, aes(x = fitness1, y = fitness2)) +
    geom_point(color = "#75C2F6", alpha = point_alpha, size = 1.5) +
    coord_cartesian(xlim = c(-1.5, 1.0), ylim = c(-1.5, 0.7)) +   # Fixed coordinate range
    labs(
      title = main_title,
      x = x_lab,
      y = y_lab
    ) +
    annotate("text", 
             x = -1.4, y = 0.45,
             label = paste0("r = ", r_value, "\np = ", 
                            ifelse(p_value < 0.0001, "< 0.0001", p_value)),
             hjust = 0, vjust = 1, size = 3,
             color = "black") +
    theme_classic(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11),
      axis.text = element_text(size = 10, colour = "black"),
      # 添加这行：旋转X轴刻度标签90度
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title = element_text(size = 10),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  # Show plot
  print(p)
  
  # Save plot
  if (!is.null(output_file)) {
    ggsave(filename = output_file,
           plot = p,
           device = cairo_pdf,
           width = plot_width,
           height = plot_height,
           units = "in",
           dpi = 300)
    cat(paste("Plot saved to:", output_file, "\n"))
  }
  
  return(list(data = data, plot = p))
}




# Wild-type amino acid sequence
wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

### Abundance - Overall only (block1 removed)
result_abundance <- compare_fitness_libraries_singlemut_overall_no_block1(
  # nicking library data
  lib1_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_abundance_2_fitness_replicates_fullseq.RData", 
  lib1_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
  
  # synthetic library data
  lib2_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  lib2_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
  
  wt_aa = wt_aa,
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figures1/Comparison_of_fitness_data_Abundance_Overall_no_block1.pdf",
  x_lab = "Abundance nicking library fitness",
  y_lab = "Abundance synthetic library fitness",
  main_title = "Comparison of fitness data between synthetic library and nicking library\nsingle mutation",
  point_alpha = 0.3,
  plot_width = 5,
  plot_height = 5
)

### RAF1 - Overall only (block1 removed)
result_raf1 <- compare_fitness_libraries_singlemut_overall_no_block1(
  # nicking library data
  lib1_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData", 
  lib1_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData",
  
  # synthetic library data
  lib2_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  lib2_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/RAF_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  
  wt_aa = wt_aa,
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figures1/Comparison_of_fitness_data_RAF1_Overall_no_block1.pdf",
  x_lab = "RAF1 nicking library fitness",
  y_lab = "RAF1 synthetic library fitness",
  main_title = "Comparison of fitness data between synthetic library and nicking library\nsingle mutation",
  point_alpha = 0.3,
  plot_width = 5,
  plot_height = 5
)

### K55 - Overall only (block1 removed)
result_K55 <- compare_fitness_libraries_singlemut_overall_no_block1(
  # nicking library data
  lib1_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K55_2_fitness_replicates_fullseq.RData", 
  lib1_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K55_3_fitness_replicates_fullseq.RData",
  
  # synthetic library data
  lib2_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K55_block2_Q20_rbg2_filter1_fitness_replicates.RData",
  lib2_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K55_block3_Q20_rbg_1_filter1_fitness_replicates.RData",
  
  wt_aa = wt_aa,
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figures1/Comparison_of_fitness_data_K55_Overall_no_block1.pdf",
  x_lab = "K55 nicking library fitness",
  y_lab = "K55 synthetic library fitness",
  main_title = "Comparison of fitness data between synthetic library and nicking library\nsingle mutation",
  point_alpha = 0.3,
  plot_width = 5,
  plot_height = 5
)

### K27 - Overall only (block1 removed)
result_K27 <- compare_fitness_libraries_singlemut_overall_no_block1(
  # nicking library data
  lib1_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K27_2_fitness_replicates_fullseq.RData", 
  lib1_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K27_3_fitness_replicates_fullseq.RData",
  
  # synthetic library data
  lib2_block2 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K27_block2_Q20_rbg3_filter1_fitness_replicates.RData",
  lib2_block3 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K27_block3_Q20_rbg3_filter1_fitness_replicates.RData",
  
  wt_aa = wt_aa,
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figures1/Comparison_of_fitness_data_K27_Overall_no_block1.pdf",
  x_lab = "K27 nicking library fitness",
  y_lab = "K27 synthetic library fitness",
  main_title = "Comparison of fitness data between synthetic library and nicking library\nsingle mutation",
  point_alpha = 0.3,
  plot_width = 5,
  plot_height = 5
)
