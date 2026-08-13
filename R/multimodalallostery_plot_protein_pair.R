#' Plot Protein Pair.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param analysis_result Value supplied for `analysis_result`.
#' @param color_map Named colors for classification groups.
#' @param legend_order Ordered classification labels used in the plot.
#' @param point_size Value supplied for `point_size`.
#' @param alpha Value supplied for `alpha`.
#' @param base_size Value supplied for `base_size`.
#' @param xlim Value supplied for `xlim`.
#' @param ylim Value supplied for `ylim`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_protein_pair <- function(analysis_result, color_map, legend_order, point_size = 2.5, alpha = 0.7,
                              base_size = 15, xlim = c(-1.5, 3), ylim = c(-1.5, 3)) {
    merged_data <- analysis_result$data
    protein_x <- analysis_result$names[1]
    protein_y <- analysis_result$names[2]
    threshold_x <- analysis_result$thresholds[1]
    threshold_y <- analysis_result$thresholds[2]
    ddG_x_col <- paste0("ddG_", protein_x)
    ddG_y_col <- paste0("ddG_", protein_y)
    cor_test <- cor.test(merged_data[[ddG_x_col]], merged_data[[ddG_y_col]])
    r_value <- round(cor_test$estimate, 3)
    p_value <- cor_test$p.value
    sig_label <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "")))
    merged_data$final_classification <- factor(merged_data$final_classification, levels = legend_order)
    anticorrelated_nbp <- merged_data[final_classification %in% c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y") & 
        is_NBP == TRUE]
    anticorrelated_non_nbp <- merged_data[final_classification %in% c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y") & 
        is_NBP == FALSE]
    other_significant <- merged_data[final_classification %in% c("Both promoting", "Both disrupting", "Allosteric only in X", 
        "Allosteric only in Y")]
    not_significant <- merged_data[final_classification == "Not significant (FDR >= 0.05)"]
    p <- ggplot() + theme_classic(base_size = base_size) + geom_vline(xintercept = c(-threshold_x, threshold_x), linetype = "dashed", 
        color = "grey50", linewidth = 0.5) + geom_hline(yintercept = c(-threshold_y, threshold_y), linetype = "dashed", color = "grey50", 
        linewidth = 0.5) + geom_point(data = not_significant, aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"), 
        size = point_size, alpha = alpha * 0.3, stroke = 0.3) + geom_point(data = other_significant, aes_string(x = ddG_x_col, 
        y = ddG_y_col, color = "final_classification"), size = point_size, alpha = alpha, stroke = 0.3) + geom_point(data = anticorrelated_non_nbp, 
        aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"), size = point_size, alpha = alpha, shape = 16, 
        stroke = 0.3) + geom_point(data = anticorrelated_nbp, aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"), 
        size = point_size + 0.5, alpha = alpha, shape = 17, stroke = 0.3) + scale_color_manual(values = color_map, breaks = legend_order, 
        drop = FALSE) + annotate("text", x = -Inf, y = Inf, label = paste0("R = ", r_value, sig_label), hjust = -0.2, vjust = 1.5, 
        size = base_size/3) + xlab(paste0("Binding \u0394\u0394G (", protein_x, ") (kcal/mol)")) + ylab(paste0("Binding \u0394\u0394G (", 
        protein_y, ") (kcal/mol)")) + theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", 
        color = NA), legend.position = "bottom", legend.text = element_text(size = base_size - 2), legend.title = element_blank(), 
        legend.key.size = unit(0.4, "cm"), legend.spacing.y = unit(0.1, "cm"), legend.margin = margin(t = 5, b = 5), axis.text = element_text(size = base_size - 
            2), axis.title = element_text(size = base_size), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), plot.margin = margin(10, 10, 10, 10)) + 
        coord_cartesian(xlim = xlim, ylim = ylim) + guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 3)))
    return(p)
}

