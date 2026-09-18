#' plot site scatter
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_site_scatter()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param full_analysis_result Input or option used by this workflow; see Usage for the exact interface.
#' @param target_position Input or option used by this workflow; see Usage for the exact interface.
#' @param point_size Input or option used by this workflow; see Usage for the exact interface.
#' @param alpha Input or option used by this workflow; see Usage for the exact interface.
#' @param base_size Input or option used by this workflow; see Usage for the exact interface.
#' @param show_labels Input or option used by this workflow; see Usage for the exact interface.
#' @param label_all Input or option used by this workflow; see Usage for the exact interface.
#' @param xlim Input or option used by this workflow; see Usage for the exact interface.
#' @param ylim Input or option used by this workflow; see Usage for the exact interface.
#' @param show_WT Input or option used by this workflow; see Usage for the exact interface.
#' @param WT_point_size Input or option used by this workflow; see Usage for the exact interface.
#' @param WT_color Input or option used by this workflow; see Usage for the exact interface.
#' @param WT_shape Input or option used by this workflow; see Usage for the exact interface.
#' @param WT_label Input or option used by this workflow; see Usage for the exact interface.
#' @param fixed_threshold Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_site_scatter <- function (full_analysis_result, target_position, point_size = 4, alpha = 0.80000000000000004, base_size = 14, show_labels = TRUE, 
    label_all = TRUE, xlim = NULL, ylim = NULL, show_WT = TRUE, WT_point_size = 3, WT_color = "black", WT_shape = 19, WT_label = "WT", 
    fixed_threshold = 0.40000000000000002) 
{
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
    cor_test <- stats::cor.test(cor_data[[paste0("ddG_", assay_x)]], cor_data[[paste0("ddG_", assay_y)]])
    r_value <- round(cor_test$estimate, 3)
    p_value <- cor_test$p.value
    cor_data_with_WT <- site_data
    cor_test_with_WT <- stats::cor.test(cor_data_with_WT[[paste0("ddG_", assay_x)]], cor_data_with_WT[[paste0("ddG_", assay_y)]])
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
        cat("\n\u6807\u8bb0\u6240\u6709\u7a81\u53d8\uff0c\u5171", nrow(data_to_label), "\u4e2a\n")
    }
    else {
        anticorrelated_groups <- c("Anticorrelated_NBP", "Anticorrelated_SwitchI", "Anticorrelated_SwitchII", "Anticorrelated_Other")
        data_to_label <- site_data[plot_group %in% anticorrelated_groups]
        cat("\n\u53ea\u6807\u8bb0anticorrelated\u7a81\u53d8\uff0c\u5171", nrow(data_to_label), "\u4e2a\n")
    }
    wt_data <- site_data[mt_codon == "WT"]
    mut_data <- site_data[mt_codon != "WT"]
    p <- ggplot2::ggplot() + ggplot2::theme_classic(base_size = base_size) + ggplot2::geom_vline(xintercept = c(-fixed_threshold, 
        fixed_threshold), linetype = "dashed", color = "grey60", linewidth = 0.80000000000000004) + ggplot2::geom_hline(yintercept = c(-fixed_threshold, 
        fixed_threshold), linetype = "dashed", color = "grey60", linewidth = 0.80000000000000004) + ggplot2::geom_point(data = mut_data[final_classification == 
        "Not significant (FDR >= 0.05)"], ggplot2::aes(x = .data[[paste0("ddG_", assay_x)]], y = .data[[paste0("ddG_", assay_y)]], 
        color = final_classification, shape = plot_shape), size = point_size * 0.69999999999999996, alpha = alpha * 0.29999999999999999, 
        stroke = 0.29999999999999999) + ggplot2::geom_point(data = mut_data[final_classification != "Not significant (FDR >= 0.05)"], 
        ggplot2::aes(x = .data[[paste0("ddG_", assay_x)]], y = .data[[paste0("ddG_", assay_y)]], color = final_classification, 
            shape = plot_shape), size = point_size, alpha = alpha, stroke = 0.80000000000000004)
    if (show_WT && nrow(wt_data) > 0) {
        p <- p + ggplot2::geom_point(data = wt_data, ggplot2::aes(x = .data[[paste0("ddG_", assay_x)]], y = .data[[paste0("ddG_", 
            assay_y)]], shape = plot_shape), color = WT_color, size = WT_point_size, stroke = 1.2)
        if (show_labels) {
            p <- p + ggrepel::geom_text_repel(data = wt_data, ggplot2::aes(x = .data[[paste0("ddG_", assay_x)]], y = .data[[paste0("ddG_", 
                assay_y)]], label = WT_label), color = WT_color, size = 4, fontface = "bold", box.padding = 0.29999999999999999, 
                point.padding = 0.20000000000000001, segment.color = "grey30", segment.size = 0.5, segment.alpha = 0.80000000000000004, 
                min.segment.length = 0, force = 1, force_pull = 0.5, seed = 123, show.legend = FALSE)
        }
    }
    p <- p + ggplot2::scale_color_manual(values = color_values, breaks = names(color_values), drop = FALSE) + ggplot2::scale_shape_manual(values = shape_values, 
        breaks = names(shape_values), drop = FALSE) + ggplot2::annotate("text", x = -Inf, y = Inf, label = paste0("R = ", 
        r_value_with_WT, ifelse(p_value_with_WT < 0.001, "***", ifelse(p_value_with_WT < 0.01, "**", ifelse(p_value_with_WT < 
            0.050000000000000003, "*", " ns")))), hjust = -0.10000000000000001, vjust = 1.5, size = base_size/3.5) + ggplot2::labs(x = bquote(Binding ~ 
        Delta * Delta * G ~ "(" * .(assay_x) * ") (kcal/mol)"), y = bquote(Binding ~ Delta * Delta * G ~ "(" * .(assay_y) * 
        ") (kcal/mol)"), title = paste0(assay_x, " vs ", assay_y, " - Position ", target_position)) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", 
        color = NA), plot.background = ggplot2::element_rect(fill = "white", color = NA), legend.position = "bottom", legend.text = ggplot2::element_text(size = base_size - 
        2), legend.title = ggplot2::element_blank(), legend.key.size = grid::unit(0.40000000000000002, "cm"), legend.spacing.y = grid::unit(0.10000000000000001, 
        "cm"), legend.margin = ggplot2::margin(t = 5, b = 5), axis.text = ggplot2::element_text(size = base_size - 2), axis.text.x = ggplot2::element_text(angle = 90, 
        hjust = 1, vjust = 0.5), axis.title = ggplot2::element_text(size = base_size), panel.border = ggplot2::element_rect(color = "black", 
        fill = NA, linewidth = 0.80000000000000004), plot.margin = ggplot2::margin(10, 10, 10, 10), plot.title = ggplot2::element_text(hjust = 0.5, 
        size = base_size + 2), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = base_size - 1, color = "grey40"))
    if (!is.null(xlim) && !is.null(ylim)) {
        p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, clip = "off")
    }
    else {
        x_range <- range(site_data[[paste0("ddG_", assay_x)]], na.rm = TRUE)
        y_range <- range(site_data[[paste0("ddG_", assay_y)]], na.rm = TRUE)
        x_pad <- diff(x_range) * 0.14999999999999999
        y_pad <- diff(y_range) * 0.14999999999999999
        p <- p + ggplot2::coord_cartesian(xlim = c(x_range[1] - x_pad, x_range[2] + x_pad), ylim = c(y_range[1] - y_pad, 
            y_range[2] + y_pad), clip = "off")
    }
    if (show_labels && nrow(data_to_label) > 0) {
        p <- p + ggrepel::geom_text_repel(data = data_to_label, ggplot2::aes(x = .data[[paste0("ddG_", assay_x)]], y = .data[[paste0("ddG_", 
            assay_y)]], label = mt, color = final_classification), size = 3.5, box.padding = 0.29999999999999999, point.padding = 0.20000000000000001, 
            segment.color = "grey50", segment.size = 0.29999999999999999, segment.alpha = 0.59999999999999998, min.segment.length = 0, 
            max.overlaps = Inf, force = 1, force_pull = 0.5, seed = 123, show.legend = FALSE)
    }
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 3)), shape = "none")
    return(p)
}

