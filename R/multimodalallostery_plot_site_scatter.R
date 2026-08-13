#' Plot Site Scatter.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param full_analysis_result Value supplied for `full_analysis_result`.
#' @param target_position Value supplied for `target_position`.
#' @param point_size Value supplied for `point_size`.
#' @param alpha Value supplied for `alpha`.
#' @param base_size Value supplied for `base_size`.
#' @param show_labels Value supplied for `show_labels`.
#' @param label_all Value supplied for `label_all`.
#' @param xlim Value supplied for `xlim`.
#' @param ylim Value supplied for `ylim`.
#' @param show_WT Value supplied for `show_WT`.
#' @param WT_point_size Value supplied for `WT_point_size`.
#' @param WT_color Value supplied for `WT_color`.
#' @param WT_shape Value supplied for `WT_shape`.
#' @param WT_label Value supplied for `WT_label`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_site_scatter <- function(full_analysis_result, target_position, point_size = 4, alpha = 0.8, base_size = 14, show_labels = TRUE, 
    label_all = TRUE, xlim = NULL, ylim = NULL, show_WT = TRUE, WT_point_size = 3, WT_color = "black", WT_shape = 19, WT_label = "WT") {
    full_data <- full_analysis_result$data
    site_data <- full_data[Pos_real == target_position]
    if (nrow(site_data) == 0) {
        stop(paste("Position", target_position, "not found in the data"))
    }
    assay_x <- full_analysis_result$assays[1]
    assay_y <- full_analysis_result$assays[2]
    threshold_x <- full_analysis_result$thresholds[assay_x]
    threshold_y <- full_analysis_result$thresholds[assay_y]
    cor_data <- site_data[mt_codon != "WT"]
    cor_test <- cor.test(cor_data[[paste0("ddG_", assay_x)]], cor_data[[paste0("ddG_", assay_y)]])
    r_value <- round(cor_test$estimate, 3)
    p_value <- cor_test$p.value
    cor_data_with_WT <- site_data
    cor_test_with_WT <- cor.test(cor_data_with_WT[[paste0("ddG_", assay_x)]], cor_data_with_WT[[paste0("ddG_", assay_y)]])
    r_value_with_WT <- round(cor_test_with_WT$estimate, 3)
    p_value_with_WT <- cor_test_with_WT$p.value
    site_data[, `:=`(plot_shape, "regular")]
    site_data[final_classification == "Not significant (FDR >= 0.05)", `:=`(plot_shape, "Not significant")]
    site_data[final_classification == "WT", `:=`(plot_shape, "WT")]
    other_types <- c("Both promoting", "Both disrupting", "Allosteric only in X", "Allosteric only in Y")
    site_data[final_classification %in% other_types, `:=`(plot_shape, "other_significant")]
    anticorrelated_types <- c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y")
    site_data[final_classification %in% anticorrelated_types & plot_group == "Anticorrelated_NBP", `:=`(plot_shape, "NBP")]
    site_data[final_classification %in% anticorrelated_types & plot_group == "Anticorrelated_SwitchI", `:=`(plot_shape, "SwitchI")]
    site_data[final_classification %in% anticorrelated_types & plot_group == "Anticorrelated_SwitchII", `:=`(plot_shape, 
        "SwitchII")]
    site_data[final_classification %in% anticorrelated_types & plot_group == "Anticorrelated_Other", `:=`(plot_shape, "Other_anticorrelated")]
    shape_values <- c(regular = 16, `Not significant` = 16, other_significant = 16, NBP = 17, SwitchI = 18, SwitchII = 15, 
        Other_anticorrelated = 8, WT = WT_shape)
    color_values <- c(`Not significant (FDR >= 0.05)` = "grey80", `Other (neutral in both)` = "grey80", `Both promoting` = "#FFB0A5", 
        `Both disrupting` = "#F4270C", `Promoting in X / Disrupting in Y` = "#F4AD0C", `Disrupting in X / Promoting in Y` = "#F1DD10", 
        `Allosteric only in X` = "#1B38A6", `Allosteric only in Y` = "#75C2F6")
    if (label_all) {
        data_to_label <- site_data[mt_codon != "WT"]
        cat("\nMark all mutations, total", nrow(data_to_label), "\u4E2A\n")
    }
    else {
        anticorrelated_groups <- c("Anticorrelated_NBP", "Anticorrelated_SwitchI", "Anticorrelated_SwitchII", "Anticorrelated_Other")
        data_to_label <- site_data[plot_group %in% anticorrelated_groups]
        cat("\nOnly anticorrelated mutations are marked, total", nrow(data_to_label), "\u4E2A\n")
    }
    wt_data <- site_data[mt_codon == "WT"]
    mut_data <- site_data[mt_codon != "WT"]
    p <- ggplot() + theme_classic(base_size = base_size) + geom_vline(xintercept = c(-threshold_x, threshold_x), linetype = "dashed", 
        color = "grey60", linewidth = 0.8) + geom_hline(yintercept = c(-threshold_y, threshold_y), linetype = "dashed", color = "grey60", 
        linewidth = 0.8) + geom_point(data = mut_data[final_classification == "Not significant (FDR >= 0.05)"], aes(x = .data[[paste0("ddG_", 
        assay_x)]], y = .data[[paste0("ddG_", assay_y)]], color = final_classification, shape = plot_shape), size = point_size * 
        0.7, alpha = alpha * 0.3, stroke = 0.3) + geom_point(data = mut_data[final_classification != "Not significant (FDR >= 0.05)"], 
        aes(x = .data[[paste0("ddG_", assay_x)]], y = .data[[paste0("ddG_", assay_y)]], color = final_classification, shape = plot_shape), 
        size = point_size, alpha = alpha, stroke = 0.8)
    if (show_WT && nrow(wt_data) > 0) {
        p <- p + geom_point(data = wt_data, aes(x = .data[[paste0("ddG_", assay_x)]], y = .data[[paste0("ddG_", assay_y)]], 
            shape = plot_shape), color = WT_color, size = WT_point_size, stroke = 1.2)
        if (show_labels) {
            p <- p + geom_text_repel(data = wt_data, aes(x = .data[[paste0("ddG_", assay_x)]], y = .data[[paste0("ddG_", 
                assay_y)]], label = WT_label), color = WT_color, size = 4, fontface = "bold", box.padding = 0.3, point.padding = 0.2, 
                segment.color = "grey30", segment.size = 0.5, segment.alpha = 0.8, min.segment.length = 0, force = 1, force_pull = 0.5, 
                seed = 123, show.legend = FALSE)
        }
    }
    p <- p + scale_color_manual(values = color_values, breaks = names(color_values), drop = FALSE) + scale_shape_manual(values = shape_values, 
        breaks = names(shape_values), drop = FALSE) + annotate("text", x = -Inf, y = Inf, label = paste0("R = ", r_value_with_WT, 
        ifelse(p_value_with_WT < 0.001, "***", ifelse(p_value_with_WT < 0.01, "**", ifelse(p_value_with_WT < 0.05, "*", " ns")))), 
        hjust = -0.1, vjust = 1.5, size = base_size/3.5) + labs(x = bquote(Binding ~ Delta * Delta * G ~ "(" * .(assay_x) * 
        ") (kcal/mol)"), y = bquote(Binding ~ Delta * Delta * G ~ "(" * .(assay_y) * ") (kcal/mol)"), title = paste0(assay_x, 
        " vs ", assay_y, " - Position ", target_position)) + theme(panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA), legend.position = "bottom", legend.text = element_text(size = base_size - 
            2), legend.title = element_blank(), legend.key.size = unit(0.4, "cm"), legend.spacing.y = unit(0.1, "cm"), legend.margin = margin(t = 5, 
            b = 5), axis.text = element_text(size = base_size - 2), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title = element_text(size = base_size), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), 
        plot.margin = margin(10, 10, 10, 10), plot.title = element_text(hjust = 0.5, size = base_size + 2), plot.subtitle = element_text(hjust = 0.5, 
            size = base_size - 1, color = "grey40"))
    if (!is.null(xlim) && !is.null(ylim)) {
        p <- p + coord_cartesian(xlim = xlim, ylim = ylim, clip = "off")
    }
    else {
        x_range <- range(site_data[[paste0("ddG_", assay_x)]], na.rm = TRUE)
        y_range <- range(site_data[[paste0("ddG_", assay_y)]], na.rm = TRUE)
        x_pad <- diff(x_range) * 0.15
        y_pad <- diff(y_range) * 0.15
        p <- p + coord_cartesian(xlim = c(x_range[1] - x_pad, x_range[2] + x_pad), ylim = c(y_range[1] - y_pad, y_range[2] + 
            y_pad), clip = "off")
    }
    if (show_labels && nrow(data_to_label) > 0) {
        p <- p + geom_text_repel(data = data_to_label, aes(x = .data[[paste0("ddG_", assay_x)]], y = .data[[paste0("ddG_", 
            assay_y)]], label = mt, color = final_classification), size = 3.5, box.padding = 0.3, point.padding = 0.2, segment.color = "grey50", 
            segment.size = 0.3, segment.alpha = 0.6, min.segment.length = 0, max.overlaps = Inf, force = 1, force_pull = 0.5, 
            seed = 123, show.legend = FALSE)
    }
    p <- p + guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 3)), shape = "none")
    return(p)
}

