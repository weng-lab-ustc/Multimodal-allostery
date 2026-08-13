#' Plot Dd g Beeswarm.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param ddG_file Value supplied for `ddG_file`.
#' @param assay_sele Selected assay identifier.
#' @param residues Value supplied for `residues`.
#' @param output_file Optional output file path.
#' @param width Value supplied for `width`.
#' @param height Value supplied for `height`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_dd_g_beeswarm <- function(ddG_file, assay_sele, residues, output_file, width = 14, height = 8) {
    ddG <- krasddpcams__read_ddG(ddG_file, assay_sele)
    ddG <- ddG[, c(1:3, 23, 26, 27)]
    ddG$Chemotype <- sapply(ddG$mt_codon, multimodalallostery_define_chemotype)
    ddG$sites <- paste(ddG$wt_codon, ddG$Pos_real, sep = "")
    df_plot <- ddG[ddG$sites %in% residues, ]
    chemotype.cols <- c(Aromatic = wes_palette("Darjeeling2", 6, type = "continuous")[1], Aliphatic = wes_palette("Darjeeling2", 
        6, type = "continuous")[2], `Polar uncharged` = wes_palette("Darjeeling2", 6, type = "continuous")[3], Positive = wes_palette("Darjeeling2", 
        6, type = "continuous")[4], Negative = wes_palette("Darjeeling2", 6, type = "continuous")[5], Special = wes_palette("Darjeeling2", 
        6, type = "continuous")[6])
    cairo_pdf(output_file, width = width, height = height)
    par(mar = c(8, 6, 2, 2))
    residues_ordered <- rev(residues)
    beeswarm.out <- split(df_plot$`mean_kcal/mol`, df_plot$sites)
    boxplot(beeswarm.out, col = "white", border = "white", outline = F, horizontal = F, ylim = c(-0.5, 2.5), frame = F, xaxt = "n", 
        yaxt = "n", xlab = "", ylab = "")
    abline(h = 0, col = "black", lty = "dotted", lwd = 2)
    site_means <- tapply(df_plot$`mean_kcal/mol`, df_plot$sites, mean, na.rm = TRUE)
    site_means <- site_means[residues_ordered]
    for (i in seq_along(residues_ordered)) {
        segments(x0 = i - 0.3, x1 = i + 0.3, y0 = site_means[i], y1 = site_means[i], col = alpha("grey40", 0.7), lwd = 4)
    }
    beeswarm_result <- beeswarm(df_plot$`mean_kcal/mol` ~ factor(df_plot$sites, levels = residues_ordered), method = "swarm", 
        do.plot = FALSE)
    plot_data <- data.frame(x = beeswarm_result$x, y = beeswarm_result$y, site = beeswarm_result$x.orig, value = beeswarm_result$y.orig)
    plot_data <- merge(plot_data, df_plot[, c("sites", "mean_kcal/mol", "mt_codon", "Chemotype")], by.x = c("site", "value"), 
        by.y = c("sites", "mean_kcal/mol"), all.x = TRUE)
    plot_data$color <- chemotype.cols[plot_data$Chemotype]
    text(x = plot_data$x, y = plot_data$y, labels = plot_data$mt_codon, col = plot_data$color, cex = 2)
    axis(1, at = 1:length(residues_ordered), labels = residues_ordered, las = 2, cex.axis = 1.2)
    mtext(side = 2, line = 3, cex = 1.5, expression(Delta * Delta * "G (kcal/mol)"))
    axis(2, las = 2, cex.axis = 1.2)
    legend("bottom", horiz = TRUE, legend = names(chemotype.cols), pch = 21, pt.bg = chemotype.cols, bty = "n", cex = 1.2, 
        inset = c(0, -0.2), xpd = TRUE)
    dev.off()
    cat("Plot saved to:", output_file, "\n")
    return(invisible())
}

