library(krasddpcams)
library(data.table)
library(dplyr)
library(beeswarm)
library(wesanderson)
library(ggplot2)
library(scales)

# Function to define the chemical type of an amino acid
define_chemotype <- function(aa) {
  aromatic <- c("F","W","Y")
  aliphatic <- c("A","V","I","L","M")
  polar_uncharged <- c("S","T","N","Q","C")
  positive <- c("K","R","H")
  negative <- c("D","E")
  special <- c("G","P")
  
  if (aa %in% aromatic) return("Aromatic")
  if (aa %in% aliphatic) return("Aliphatic")
  if (aa %in% polar_uncharged) return("Polar uncharged")
  if (aa %in% positive) return("Positive")
  if (aa %in% negative) return("Negative")
  if (aa %in% special) return("Special")
  return(NA)
}

# Main plotting function
plot_ddG_beeswarm <- function(ddG_file, assay_sele, residues, output_file, 
                              width = 14, height = 8) {
  
  ## ===== 1. Reading and processing data =====
  ddG <- krasddpcams__read_ddG(ddG_file, assay_sele)
  ddG <- ddG[, c(1:3, 23, 26, 27)]
  
  # Add chemical type
  ddG$Chemotype <- sapply(ddG$mt_codon, define_chemotype)
  ddG$sites <- paste(ddG$wt_codon, ddG$Pos_real, sep = "")
  
  ## ===== 2. Screening interface residues =====
  df_plot <- ddG[ddG$sites %in% residues, ]
  
  # Color by chemotype
  chemotype.cols <- c(
    "Aromatic"       = wes_palette("Darjeeling2", 6, type = "continuous")[1],
    "Aliphatic"      = wes_palette("Darjeeling2", 6, type = "continuous")[2],
    "Polar uncharged"= wes_palette("Darjeeling2", 6, type = "continuous")[3],
    "Positive"       = wes_palette("Darjeeling2", 6, type = "continuous")[4],
    "Negative"       = wes_palette("Darjeeling2", 6, type = "continuous")[5],
    "Special"        = wes_palette("Darjeeling2", 6, type = "continuous")[6]
  )
  
  ## ===== 3. Plotting =====
  #
  cairo_pdf(output_file, width = width, height = height)
  par(mar = c(8, 6, 2, 2))
  #
  residues_ordered <- rev(residues)
  # 
  beeswarm.out <- split(df_plot$`mean_kcal/mol`, df_plot$sites)
  #
  boxplot(beeswarm.out, col = "white", border = "white", outline = F, horizontal = F,
          ylim = c(-0.5, 2.5),
          frame = F, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
  # 
  abline(h = 0, col = "black", lty = "dotted", lwd = 2)
  # 
  site_means <- tapply(df_plot$`mean_kcal/mol`, df_plot$sites, mean, na.rm = TRUE)
  site_means <- site_means[residues_ordered]  
  for (i in seq_along(residues_ordered)) {
    segments(x0 = i - 0.3, x1 = i + 0.3, 
             y0 = site_means[i],
             y1 = site_means[i], 
             col = alpha("grey40", 0.7), lwd = 4)
  }
  beeswarm_result <- beeswarm(df_plot$`mean_kcal/mol` ~ factor(df_plot$sites, levels = residues_ordered),
                              method = "swarm", 
                              do.plot = FALSE)
  plot_data <- data.frame(
    x = beeswarm_result$x,
    y = beeswarm_result$y,
    site = beeswarm_result$x.orig,
    value = beeswarm_result$y.orig
  )
  
  plot_data <- merge(plot_data, df_plot[, c("sites", "mean_kcal/mol", "mt_codon", "Chemotype")],
                     by.x = c("site", "value"), 
                     by.y = c("sites", "mean_kcal/mol"),
                     all.x = TRUE)
  plot_data$color <- chemotype.cols[plot_data$Chemotype]
  text(x = plot_data$x, 
       y = plot_data$y, 
       labels = plot_data$mt_codon, 
       col = plot_data$color, 
       cex = 2)
  # 
  axis(1, at = 1:length(residues_ordered), labels = residues_ordered, las = 2, cex.axis = 1.2)
  # 
  mtext(side = 2, line = 3, cex = 1.5,
        expression(Delta * Delta * "G (kcal/mol)"))
  # 
  axis(2, las = 2, cex.axis = 1.2)
  # 
  legend("bottom", horiz = TRUE,
         legend = names(chemotype.cols),
         pch = 21, pt.bg = chemotype.cols,
         bty = "n", cex = 1.2, inset = c(0, -0.2), xpd = TRUE)
  
  dev.off()
  
  cat("Plot saved to:", output_file, "\n")
  return(invisible())
}

##### Usage
residues_list <- c("K88", "E91", "T87", "Q129", "F90", "L133", "H94", 
                   "Y137", "H95", "R68", "S136", "Q99", "R102", "K101", 
                   "E107", "E98")

plot_ddG_beeswarm(
  ddG_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  assay_sele = "K13",
  residues = residues_list,
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figure 2/Figure_sites_mut_ddG_beeswarm_K13.pdf"
)



plot_ddG_beeswarm(
  ddG_file = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt", 
  assay_sele = "K19",
  residues = residues_list,
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figure 2/Figure_sites_mut_ddG_beeswarm_K19.pdf"
)
