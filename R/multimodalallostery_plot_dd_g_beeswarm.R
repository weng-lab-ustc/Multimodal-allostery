#' plot dd g beeswarm
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_dd_g_beeswarm()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param ddG_file Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_sele Input or option used by this workflow; see Usage for the exact interface.
#' @param residues Input or option used by this workflow; see Usage for the exact interface.
#' @param output_file Input or option used by this workflow; see Usage for the exact interface.
#' @param width Input or option used by this workflow; see Usage for the exact interface.
#' @param height Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_dd_g_beeswarm <- function (ddG_file, assay_sele, residues, output_file, width = 14, height = 8) 
{
    ddG <- krasddpcams::krasddpcams__read_ddG(ddG_file, assay_sele)
    ddG <- ddG[, c(1:3, 23, 26, 27)]
    ddG$Chemotype <- sapply(ddG$mt_codon, multimodalallostery_define_chemotype)
    ddG$sites <- paste(ddG$wt_codon, ddG$Pos_real, sep = "")
    df_plot <- ddG[ddG$sites %in% residues, ]
    chemotype.cols <- c(Aromatic = wesanderson::wes_palette("Darjeeling2", 6, type = "continuous")[1], Aliphatic = wesanderson::wes_palette("Darjeeling2", 
        6, type = "continuous")[2], `Polar uncharged` = wesanderson::wes_palette("Darjeeling2", 6, type = "continuous")[3], 
        Positive = wesanderson::wes_palette("Darjeeling2", 6, type = "continuous")[4], Negative = wesanderson::wes_palette("Darjeeling2", 
            6, type = "continuous")[5], Special = wesanderson::wes_palette("Darjeeling2", 6, type = "continuous")[6])
    grDevices::cairo_pdf(output_file, width = width, height = height)
    graphics::par(mar = c(8, 6, 2, 2))
    residues_ordered <- rev(residues)
    beeswarm.out <- split(df_plot$`mean_kcal/mol`, df_plot$sites)
    graphics::boxplot(beeswarm.out, col = "white", border = "white", outline = F, horizontal = F, ylim = c(-0.5, 2.5), frame = F, 
        xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    graphics::abline(h = 0, col = "black", lty = "dotted", lwd = 2)
    site_means <- tapply(df_plot$`mean_kcal/mol`, df_plot$sites, mean, na.rm = TRUE)
    site_means <- site_means[residues_ordered]
    for (i in seq_along(residues_ordered)) {
        graphics::segments(x0 = i - 0.29999999999999999, x1 = i + 0.29999999999999999, y0 = site_means[i], y1 = site_means[i], 
            col = scales::alpha("grey40", 0.69999999999999996), lwd = 4)
    }
    beeswarm_result <- beeswarm::beeswarm(df_plot$`mean_kcal/mol` ~ factor(df_plot$sites, levels = residues_ordered), method = "swarm", 
        do.plot = FALSE)
    plot_data <- data.frame(x = beeswarm_result$x, y = beeswarm_result$y, site = beeswarm_result$x.orig, value = beeswarm_result$y.orig)
    plot_data <- merge(plot_data, df_plot[, c("sites", "mean_kcal/mol", "mt_codon", "Chemotype")], by.x = c("site", "value"), 
        by.y = c("sites", "mean_kcal/mol"), all.x = TRUE)
    plot_data$color <- chemotype.cols[plot_data$Chemotype]
    graphics::text(x = plot_data$x, y = plot_data$y, labels = plot_data$mt_codon, col = plot_data$color, cex = 2)
    graphics::axis(1, at = 1:length(residues_ordered), labels = residues_ordered, las = 2, cex.axis = 1.2)
    graphics::mtext(side = 2, line = 3, cex = 1.5, expression(Delta * Delta * "G (kcal/mol)"))
    graphics::axis(2, las = 2, cex.axis = 1.2)
    graphics::legend("bottom", horiz = TRUE, legend = names(chemotype.cols), pch = 21, pt.bg = chemotype.cols, bty = "n", 
        cex = 1.2, inset = c(0, -0.20000000000000001), xpd = TRUE)
    grDevices::dev.off()
    cat("Plot saved to:", output_file, "\n")
    return(invisible())
}

