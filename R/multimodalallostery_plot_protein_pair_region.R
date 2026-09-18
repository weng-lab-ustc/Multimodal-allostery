#' plot protein pair region
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_protein_pair_region()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param analysis_result Input or option used by this workflow; see Usage for the exact interface.
#' @param region_residues Input or option used by this workflow; see Usage for the exact interface.
#' @param region_name Input or option used by this workflow; see Usage for the exact interface.
#' @param point_size Input or option used by this workflow; see Usage for the exact interface.
#' @param alpha Input or option used by this workflow; see Usage for the exact interface.
#' @param base_size Input or option used by this workflow; see Usage for the exact interface.
#' @param xlim Input or option used by this workflow; see Usage for the exact interface.
#' @param ylim Input or option used by this workflow; see Usage for the exact interface.
#' @param color_map Input or option used by this workflow; see Usage for the exact interface.
#' @param legend_order Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_protein_pair_region <- function (analysis_result, region_residues, region_name, point_size = 2.5, alpha = 0.69999999999999996, base_size = 12, xlim = c(-1.5, 
    3), ylim = c(-1.5, 3), color_map = c(`Not significant (FDR >= 0.05)` = "grey90", `Other (neutral in both)` = "grey90", 
    `Both promoting` = "#FFB0A5", `Both disrupting` = "#F4270C", `Promoting in X / Disrupting in Y` = "#F4AD0C", `Disrupting in X / Promoting in Y` = "#F1DD10", 
    `Allosteric only in X` = "#1B38A6", `Allosteric only in Y` = "#75C2F6"), legend_order = c("Both promoting", "Both disrupting", 
    "Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y", "Allosteric only in X", "Allosteric only in Y", 
    "Other (neutral in both)", "Not significant (FDR >= 0.05)")) 
{
    merged_data <- analysis_result$data
    protein_x <- analysis_result$names[1]
    protein_y <- analysis_result$names[2]
    threshold_x <- analysis_result$thresholds[1]
    threshold_y <- analysis_result$thresholds[2]
    region_data <- merged_data[Pos_real %in% region_residues]
    if (nrow(region_data) == 0) {
        cat("\u8b66\u544a: \u533a\u57df", region_name, "\u6ca1\u6709\u7a81\u53d8\u6570\u636e\n")
        return(NULL)
    }
    ddG_x_col <- paste0("ddG_", protein_x)
    ddG_y_col <- paste0("ddG_", protein_y)
    if (nrow(region_data) > 2) {
        cor_test <- stats::cor.test(region_data[[ddG_x_col]], region_data[[ddG_y_col]])
        r_value <- round(cor_test$estimate, 3)
        p_value <- cor_test$p.value
        sig_label <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.050000000000000003, "*", 
            "")))
        cor_label <- paste0("R = ", r_value, sig_label)
    }
    else {
        cor_label <- "R = NA"
    }
    anticorrelated_nbp <- region_data[final_classification %in% c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y") & 
        is_NBP == TRUE]
    anticorrelated_non_nbp <- region_data[final_classification %in% c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y") & 
        is_NBP == FALSE]
    other_significant <- region_data[final_classification %in% c("Both promoting", "Both disrupting", "Allosteric only in X", 
        "Allosteric only in Y")]
    not_significant <- region_data[final_classification == "Not significant (FDR >= 0.05)"]
    p <- ggplot2::ggplot() + ggplot2::theme_classic(base_size = base_size) + ggplot2::geom_vline(xintercept = c(-threshold_x, 
        threshold_x), linetype = "dashed", color = "grey50", linewidth = 0.5) + ggplot2::geom_hline(yintercept = c(-threshold_y, 
        threshold_y), linetype = "dashed", color = "grey50", linewidth = 0.5) + ggplot2::geom_point(data = not_significant, 
        .ma_aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"), size = point_size * 0.80000000000000004, 
        alpha = alpha * 0.29999999999999999, stroke = 0.29999999999999999) + ggplot2::geom_point(data = other_significant, 
        .ma_aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"), size = point_size, alpha = alpha, stroke = 0.29999999999999999) + 
        ggplot2::geom_point(data = anticorrelated_non_nbp, .ma_aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"), 
            size = point_size, alpha = alpha, shape = 16, stroke = 0.29999999999999999) + ggplot2::geom_point(data = anticorrelated_nbp, 
        .ma_aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"), size = point_size + 0.5, alpha = alpha, 
        shape = 17, stroke = 0.29999999999999999) + ggplot2::scale_color_manual(values = color_map, breaks = legend_order, 
        drop = FALSE) + ggplot2::annotate("text", x = -Inf, y = Inf, label = cor_label, hjust = -0.10000000000000001, vjust = 1.5, 
        size = base_size/3) + ggplot2::xlab(paste0("\u0394\u0394Gb (", protein_x, ") (kcal/mol)")) + ggplot2::ylab(paste0("\u0394\u0394Gb (", 
        protein_y, ") (kcal/mol)")) + ggplot2::ggtitle(region_name) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", 
        color = NA), plot.background = ggplot2::element_rect(fill = "white", color = NA), legend.position = "none", legend.text = ggplot2::element_text(size = base_size - 
        2), legend.title = ggplot2::element_blank(), legend.key.size = grid::unit(0.29999999999999999, "cm"), axis.text = ggplot2::element_text(size = base_size - 
        2), axis.title = ggplot2::element_text(size = base_size), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, 
        vjust = 0.5), panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.80000000000000004), 
        plot.margin = ggplot2::margin(5, 5, 5, 5), plot.title = ggplot2::element_text(size = base_size + 2, hjust = 0.5)) + 
        ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
    return(p)
}

