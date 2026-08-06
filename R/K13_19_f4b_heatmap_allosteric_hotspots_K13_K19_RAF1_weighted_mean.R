library(krasddpcams)
library(data.table)
library(ggplot2)

colour_scheme <- list(
  "blue" = "#1B38A6",
  "red" = "#F4270C",
  "orange" = "#F4AD0C",
  "green" = "#09B636",
  "yellow" = "#F1DD10",
  "purple" = "#C68EFD",
  "hot pink" = "#FF0066",
  "light blue" = "#75C2F6",
  "light red" = "#FF6A56",
  "dark red" = "#A31300",
  "dark green" = "#007A20",
  "pink" = "#FFB0A5"
)

##### function (integrating three binding partners; no dashed boxes)
krasddpcams__plot_triple_ddG_heatmap <- function(
    ddG_K13, ddG_K19, ddG_RAF1, anno, wt_aa, colour_scheme,
    allosteric_sites_list = NULL,  # Named list containing the allosteric sites of K13, K19, and RAF1
    legend_limits = c(-1.3, 3)
) {
  
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  
  heatmap_tool <- data.table(
    wt_codon = rep(unlist(strsplit(wt_aa, "")), each = 20),
    Pos_real = rep(2:188, each = 20),
    mt_codon = unlist(aa_list)
  )
  
  # Load three datasets.
  ddG1 <- krasddpcams__read_ddG(ddG_K13, "K13")
  input1_heatmap <- merge(ddG1, heatmap_tool, by = c("wt_codon", "Pos_real", "mt_codon"), all = T)
  
  ddG2 <- krasddpcams__read_ddG(ddG_K19, "K19")
  input2_heatmap <- merge(ddG2, heatmap_tool, by = c("wt_codon", "Pos_real", "mt_codon"), all = T)
  
  ddG3 <- krasddpcams__read_ddG(ddG_RAF1, "RAF1")
  input3_heatmap <- merge(ddG3, heatmap_tool, by = c("wt_codon", "Pos_real", "mt_codon"), all = T)
  
  ddG <- rbind(input1_heatmap, input2_heatmap, input3_heatmap)
  
  ddG[wt_codon == mt_codon, `:=`(`mean_kcal/mol`, 0)]
  
  output <- merge(ddG, anno, by.x = "Pos_real", by.y = "Pos", all = T)
  
  # Screening for allosteric sites
  if (!is.null(allosteric_sites_list)) {
    all_sites <- unique(c(
      allosteric_sites_list$K13,
      allosteric_sites_list$K19,
      allosteric_sites_list$RAF1
    ))
    output <- output[Pos_real %in% all_sites, ]
  }
  
  # Sorted by site (N-terminus to C-terminus)
  output[, Pos_real := factor(Pos_real, levels = sort(unique(Pos_real)))]
  
  # Set the factor levels for mt_codon and assay.
  output <- within(output, mt_codon <- factor(mt_codon,
                                              levels = c("D","E","R","H","K","S","T","N","Q",
                                                         "C","G","P","A","V","I","L","M","F","W","Y")))
  
  output <- within(output, assay <- factor(assay, levels = c("K13", "K19", "RAF1")))
  
  # Create wtcodon_pos tags.
  output[, wtcodon_pos := paste0(wt_codon, Pos_real)]
  output[, Pos_num := as.numeric(as.character(Pos_real))]
  output <- output[order(Pos_num), ]
  output[, wtcodon_pos := factor(wtcodon_pos, levels = unique(wtcodon_pos))]
  
  # plot
  p <- ggplot2::ggplot() + 
    geom_tile(
      data = output,
      aes(x = assay, y = mt_codon, fill = `mean_kcal/mol`)
    ) +
    geom_text(
      data = output[Pos_num > 1 & wt_codon == mt_codon, ],
      aes(x = assay, y = mt_codon, label = "-"),
      size = 6 * 5/14
    ) +
    scale_fill_gradient2(
      limits = legend_limits,
      low = colour_scheme[["blue"]],
      mid = "gray",
      high = colour_scheme[["red"]],
      na.value = "white"
    ) +
    facet_wrap(~wtcodon_pos, nrow = 2) +
    ylab("Mutant AA") +
    xlab("Binding partners") +
    labs(fill = expression(Delta*Delta*G~"(kcal/mol)")) +
    theme_bw() +
    theme(
      text = element_text(size = 15, family = "Arial"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(colour = "white", fill = "white"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(margin = margin(0,5,0,0)),
      panel.spacing.y = unit(3, "mm"),
      panel.spacing.x = unit(1, "mm")
    ) +
    coord_fixed()
  
  return(p)
}

##### =====================
##### Data Preparation
##### =====================

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

anno <- fread("C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv")

# Allosteric site information
K13_all_allosteric_hotspots <- c(15, 145, 10, 21, 56, 130, 139, 142, 151, 155, 178)
K19_all_allosteric_hotspots <- c(15, 145, 10, 151, 157, 184)
RAF1_all_allosteric_hotspots <- c(15, 16, 17, 18, 28, 32, 34, 35, 57, 60, 145, 146, 6, 10, 20, 22, 54, 55, 58, 59, 77, 144, 163, 184)

# Merge all allosteric sites and sort them from N- to C-terminus.
all_allosteric_sites <- unique(c(K13_all_allosteric_hotspots, 
                                 K19_all_allosteric_hotspots, 
                                 RAF1_all_allosteric_hotspots))
all_allosteric_sites <- sort(all_allosteric_sites)

# Prepare a list of allosteric sites (for the function)
allosteric_list <- list(
  K13 = K13_all_allosteric_hotspots,
  K19 = K19_all_allosteric_hotspots,
  RAF1 = RAF1_all_allosteric_hotspots
)

##### =====================

p_triple <- krasddpcams__plot_triple_ddG_heatmap(
  ddG_K13 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  ddG_K19 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt",
  ddG_RAF1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
  anno = anno,
  wt_aa = wt_aa,
  colour_scheme = colour_scheme,
  allosteric_sites_list = allosteric_list,
  legend_limits = c(-1.3, 3)
)

#
print(p_triple)

##### =====================
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260406_start version/figure3/triple_allosteric_sites_heatmap2.pdf", 
       p_triple, width = 10, height = 8, device = cairo_pdf)

# 也可以保存为PNG格式方便查看
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260406_start version/figure3/triple_allosteric_sites_heatmap.png", 
       p_triple, width = 12, height = 10, dpi = 300)
