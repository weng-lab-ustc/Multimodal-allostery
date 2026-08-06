library(data.table)
library(ggplot2)
library(dplyr)

plot_binding_interface_residue_median_ddG_heatmap <- function(
    ddG_file,
    binding_sites,
    position_labels,
    title = "Binding Interface Residues - Median ΔΔGb"
) {
  ddG <- fread(ddG_file)
  ddG <- ddG[, Pos_real := Pos + 1]

  median_values <- ddG %>%
    filter(Pos_real %in% binding_sites) %>%
    group_by(Pos_real) %>%
    summarise(median_ddG = median(`mean_kcal/mol`, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      residue_label = factor(position_labels[as.character(Pos_real)], 
                             levels = position_labels),
    ) %>%
    arrange(residue_label)

  p <- ggplot(median_values, aes(x = 1, y = residue_label, fill = median_ddG)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_gradient2(
      low = "#1B38A6",
      mid = "gray",
      high = "#F4270C",
      midpoint = 0,
      limits = c(-1, 2.5),  
      name = expression(Delta * Delta * "Gb (kcal/mol)")
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Residue"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      legend.position = "right",
      panel.border = element_rect(color = "gray90", fill = NA, linewidth = 0.8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    ) +
    scale_y_discrete(limits = rev(levels(median_values$residue_label))) +
    coord_fixed(ratio = 0.4)
  
  return(p)
}

#===========================================================
## K13

#ddG_file <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt"
ddG_file <- "~/Library/CloudStorage/OneDrive-个人/文档/Data/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt"
# Define the binding interface location and the corresponding label.
binding_interface_residues <- c(63, 105, 106, 98, 107, 101, 102, 99, 136, 95, 137, 94, 133, 90, 129, 87, 91, 88)

position_labels <- c(
  "106" = "S106","105" = "D105","63" = "E63","98" = "E98", "107" = "E107", "101" = "K101", "102" = "R102",
  "99" = "Q99", "136" = "S136", "95" = "H95",
  "137" = "Y137", "94" = "H94", "133" = "L133", "90" = "F90",
  "129" = "Q129", "87" = "T87", "91" = "E91", "88" = "K88"
)

p <- plot_binding_interface_residue_median_ddG_heatmap(
  ddG_file = ddG_file,
  binding_sites = binding_interface_residues,
  position_labels = position_labels,
  title = "K13 Binding Interface Residues - Median ΔΔGb"
)


print(p)

ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result2_figure2/K13_BI_residues_median_ddG_heatmap.pdf", p , device = cairo_pdf,height = 8,width=3)




#=========================================================================
## K19

ddG_file <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt"


binding_interface_residues <- c(98, 107, 101, 102, 99, 136, 95, 137, 94, 133, 90, 129, 87, 91, 88, 68, 108)

position_labels <- c(
  "98" = "E98", "107" = "E107", "101" = "K101", "102" = "R102",
  "99" = "Q99", "136" = "S136", "95" = "H95",
  "137" = "Y137", "94" = "H94", "133" = "L133", "90" = "F90",
  "129" = "Q129", "87" = "T87", "91" = "E91", "88" = "K88", "68" = "R68","108" = "D108"
)



p <- plot_binding_interface_residue_median_ddG_heatmap(
  ddG_file = ddG_file,
  binding_sites = binding_interface_residues,
  position_labels = position_labels,
  title = "K19 Binding Interface Residues - Median ΔΔGb"
)


print(p)

ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result2_figure2/K19_BI_residues_median_ddG_heatmap.pdf", p , device = cairo_pdf,height = 7,width=3)


