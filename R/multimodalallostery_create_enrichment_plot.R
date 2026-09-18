#' create enrichment plot
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_create_enrichment_plot()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param plot_data Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_name Input or option used by this workflow; see Usage for the exact interface.
#' @param output_dir Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_create_enrichment_plot <- function (plot_data, assay_name, output_dir) 
{
    category_colors <- c(Allosteric = "#F4270C", Inhibit = "#75C2F6", Stabilize = "#F4AD0C")
    plot_data[, `:=`(x_offset, dplyr::case_when(Category == "Allosteric" ~ -0.25, Category == "Inhibit" ~ 0, Category == 
        "Stabilize" ~ 0.25))]
    plot_data[, `:=`(Proportion, Percentage/100)]
    plot_data[, `:=`(Error_prop, Error/100)]
    fixed_y_position <- max(plot_data$Proportion + plot_data$Error_prop, na.rm = TRUE) * 1.1499999999999999
    fixed_y_position <- max(fixed_y_position, max(plot_data$Proportion) * 1.2)
    fixed_y_position <- min(fixed_y_position, 0.69999999999999996)
    max_y <- 0.59999999999999998
    plot_data[, `:=`(Sig_label, ifelse(P_value < 0.050000000000000003, ifelse(P_value < 0.01, ifelse(P_value < 0.001, "***", 
        "**"), "*"), " ns"))]
    plot_data[, `:=`(Label, paste0("OR = ", round(OR, 2), Sig_label))]
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Region, y = Proportion, fill = Category)) + ggplot2::geom_col(position = ggplot2::position_dodge(0.80000000000000004), 
        width = 0.69999999999999996, color = "white", linewidth = 0.29999999999999999) + ggplot2::geom_errorbar(ggplot2::aes(ymin = Proportion - 
        Error_prop, ymax = Proportion + Error_prop), position = ggplot2::position_dodge(0.80000000000000004), width = 0.14999999999999999, 
        linewidth = 0.80000000000000004, color = "#F1DD10", alpha = 0.80000000000000004, na.rm = TRUE) + ggplot2::geom_errorbar(ggplot2::aes(ymin = Proportion - 
        Error_prop, ymax = Proportion - Error_prop), position = ggplot2::position_dodge(0.80000000000000004), width = 0.29999999999999999, 
        linewidth = 0.80000000000000004, color = "#F1DD10", alpha = 0.80000000000000004, na.rm = TRUE) + ggplot2::geom_errorbar(ggplot2::aes(ymin = Proportion + 
        Error_prop, ymax = Proportion + Error_prop), position = ggplot2::position_dodge(0.80000000000000004), width = 0.29999999999999999, 
        linewidth = 0.80000000000000004, color = "#F1DD10", alpha = 0.80000000000000004, na.rm = TRUE) + ggplot2::geom_text(ggplot2::aes(x = as.numeric(Region) + 
        x_offset, y = fixed_y_position, label = Label), size = 3.5, angle = 45, hjust = 0.5, vjust = 0, color = ifelse(plot_data$Category == 
        "Allosteric", "#F4270C", ifelse(plot_data$Category == "Inhibit", "#75C2F6", "#F4AD0C"))) + ggplot2::theme_classic(base_size = 15) + 
        ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5, margin = ggplot2::margin(b = 15)), axis.title = ggplot2::element_text(size = 14), 
            axis.text = ggplot2::element_text(size = 12, color = "black"), axis.text.x = ggplot2::element_text(size = 13), 
            axis.line = ggplot2::element_line(color = "black", linewidth = 0.5), axis.ticks = ggplot2::element_line(color = "black", 
                linewidth = 0.5), legend.title = ggplot2::element_text(size = 13), legend.text = ggplot2::element_text(size = 12), 
            legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", panel.grid = ggplot2::element_blank(), 
            plot.margin = ggplot2::margin(t = 40, r = 10, b = 10, l = 10)) + ggplot2::geom_hline(yintercept = 0, color = "black", 
        linewidth = 0.5) + ggplot2::labs(title = paste0("Enrichment of Allosteric Mutations - ", assay_name), x = NULL, y = "Proportion of all allosteric mutations", 
        fill = "Direction") + ggplot2::scale_y_continuous(limits = c(0, max_y), expand = ggplot2::expansion(mult = c(0, 0.02)), 
        breaks = seq(0, max_y, 0.10000000000000001), labels = scales::percent) + ggplot2::scale_fill_manual(values = category_colors)
    .ma_save_plot(filename = file.path(output_dir, paste0("allosteric mutations Enrichment_plot_", assay_name, ".pdf")), 
        plot = p, width = 8, height = 6, dpi = 300)
    .ma_save_plot(filename = file.path(output_dir, paste0("allosteric mutations Enrichment_plot_", assay_name, ".png")), 
        plot = p, width = 8, height = 6, dpi = 300)
    return(p)
}

