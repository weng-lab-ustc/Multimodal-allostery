#' Plot Scatter Dd Gb Dd Gf.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param ddG1 Value supplied for `ddG1`.
#' @param assay1 Value supplied for `assay1`.
#' @param ddG2 Value supplied for `ddG2`.
#' @param assay2 Value supplied for `assay2`.
#' @param anno Annotation data or annotation file path used by the analysis.
#' @param binder Value supplied for `binder`.
#' @param colour_scheme Value supplied for `colour_scheme`.
#' @param xlim Value supplied for `xlim`.
#' @param ylim Value supplied for `ylim`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_scatter_dd_gb_dd_gf <- function(ddG1, assay1, ddG2, assay2, anno, binder = c("K13", "K19", "RAF1"), colour_scheme, xlim = c(-1, 
    3.3), ylim = c(-1.7, 3)) {
    binder <- match.arg(binder)
    ddG1 <- krasddpcams__read_ddG(ddG = ddG1, assay_sele = assay1)
    ddG2 <- krasddpcams__read_ddG(ddG = ddG2, assay_sele = assay2)
    all_ddG <- rbind(ddG1, ddG2)
    all_ddG_dc <- dcast(all_ddG[!is.na(mt), ], mt + Pos_real ~ assay, value.var = "mean_kcal/mol")
    all_ddG_dc_anno <- merge(all_ddG_dc, anno, by.x = "Pos_real", by.y = "Pos")
    dist_col <- paste0("scHAmin_ligand_", binder)
    all_ddG_dc_anno[, `:=`(binding_type, "others")]
    all_ddG_dc_anno[get(dist_col) <= 5, `:=`(binding_type, "binding interface")]
    p <- ggplot() + geom_point(data = all_ddG_dc_anno[Pos_real > 1 & binding_type == "others", ], aes(x = !!sym(assay1), 
        y = !!sym(binder)), color = "grey70", alpha = 0.4, size = 1.5) + geom_point(data = all_ddG_dc_anno[Pos_real > 1 & 
        binding_type == "binding interface"], aes(x = !!sym(assay1), y = !!sym(binder)), color = "#F4270C", alpha = 0.5, 
        size = 2) + scale_x_continuous(limits = xlim, expand = c(0, 0)) + scale_y_continuous(limits = ylim, expand = c(0, 
        0)) + theme_classic() + labs(x = "Folding \u0394\u0394G (kcal mol-1)", y = paste0(binder, " Binding \u0394\u0394G (kcal mol-1)")) + 
        theme(text = element_text(size = 8), legend.position = "right", legend.text = element_text(size = 8), axis.text.x = element_text(angle = 90, 
            vjust = 0.5, hjust = 1, size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black", vjust = 0.5, 
            hjust = 0.5, margin = margin(0, 0, 0, 0, "mm")), axis.title = element_text(size = 8), legend.key.height = unit(3.1, 
            "mm"), legend.key.width = unit(3.1, "mm"), legend.key.size = unit(1, "mm"), plot.margin = margin(3, 3, 3, 3)) + 
        coord_fixed(ratio = 1, expand = TRUE, clip = "on")
    return(p)
}

