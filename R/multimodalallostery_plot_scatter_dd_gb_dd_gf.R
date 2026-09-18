#' plot scatter dd gb dd gf
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_scatter_dd_gb_dd_gf()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param ddG1 Input or option used by this workflow; see Usage for the exact interface.
#' @param assay1 Input or option used by this workflow; see Usage for the exact interface.
#' @param ddG2 Input or option used by this workflow; see Usage for the exact interface.
#' @param assay2 Input or option used by this workflow; see Usage for the exact interface.
#' @param anno Input or option used by this workflow; see Usage for the exact interface.
#' @param binder Input or option used by this workflow; see Usage for the exact interface.
#' @param colour_scheme Input or option used by this workflow; see Usage for the exact interface.
#' @param xlim Input or option used by this workflow; see Usage for the exact interface.
#' @param ylim Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_scatter_dd_gb_dd_gf <- function (ddG1, assay1, ddG2, assay2, anno, binder = c("K13", "K19", "RAF1"), colour_scheme, xlim = c(-1, 3.2999999999999998), 
    ylim = c(-1.7, 3)) 
{
    binder <- match.arg(binder)
    ddG1 <- krasddpcams::krasddpcams__read_ddG(ddG = ddG1, assay_sele = assay1)
    ddG2 <- krasddpcams::krasddpcams__read_ddG(ddG = ddG2, assay_sele = assay2)
    all_ddG <- rbind(ddG1, ddG2)
    all_ddG_dc <- data.table::dcast(all_ddG[!is.na(mt), ], mt + Pos_real ~ assay, value.var = "mean_kcal/mol")
    all_ddG_dc_anno <- merge(all_ddG_dc, anno, by.x = "Pos_real", by.y = "Pos")
    dist_col <- paste0("scHAmin_ligand_", binder)
    all_ddG_dc_anno[, `:=`(binding_type, "others")]
    all_ddG_dc_anno[get(dist_col) <= 5, `:=`(binding_type, "binding interface")]
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = all_ddG_dc_anno[Pos_real > 1 & binding_type == "others", ], ggplot2::aes(x = !!rlang::sym(assay1), 
        y = !!rlang::sym(binder)), color = "grey70", alpha = 0.40000000000000002, size = 1.5) + ggplot2::geom_point(data = all_ddG_dc_anno[Pos_real > 
        1 & binding_type == "binding interface"], ggplot2::aes(x = !!rlang::sym(assay1), y = !!rlang::sym(binder)), color = "#F4270C", 
        alpha = 0.5, size = 2) + ggplot2::scale_x_continuous(limits = xlim, expand = c(0, 0)) + ggplot2::scale_y_continuous(limits = ylim, 
        expand = c(0, 0)) + ggplot2::theme_classic() + ggplot2::labs(x = "Folding \u0394\u0394G (kcal mol-1)", y = paste0(binder, 
        " Binding \u0394\u0394G (kcal mol-1)")) + ggplot2::theme(text = ggplot2::element_text(size = 8), legend.position = "right", 
        legend.text = ggplot2::element_text(size = 8), axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, 
            size = 8, colour = "black"), axis.text.y = ggplot2::element_text(size = 8, colour = "black", vjust = 0.5, hjust = 0.5, 
            margin = ggplot2::margin(0, 0, 0, 0, "mm")), axis.title = ggplot2::element_text(size = 8), legend.key.height = grid::unit(3.1000000000000001, 
            "mm"), legend.key.width = grid::unit(3.1000000000000001, "mm"), legend.key.size = grid::unit(1, "mm"), plot.margin = ggplot2::margin(3, 
            3, 3, 3)) + ggplot2::coord_fixed(ratio = 1, expand = TRUE, clip = "on")
    return(p)
}

