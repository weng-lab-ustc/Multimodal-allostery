#' Create Enrichment Plot.
#'
#' Reusable utility logic consolidated from the selected Figure/Panel implementations.
#'
#' @param plot_data Prepared data used for plotting.
#' @param assay_name Value supplied for `assay_name`.
#' @param output_dir Value supplied for `output_dir`.
#'
#' @return The transformed result produced by the function.
#' @export
multimodalallostery_create_enrichment_plot <- function(plot_data, assay_name, output_dir) {
    category_colors <- c(Allosteric = "#F4270C", Inhibit = "#75C2F6", Stabilize = "#F4AD0C")
    plot_data[, `:=`(x_offset, dplyr::case_when(Category == "Allosteric" ~ -0.25, Category == "Inhibit" ~ 0, Category == 
        "Stabilize" ~ 0.25))]
    plot_data[, `:=`(Proportion, Percentage/100)]
    plot_data[, `:=`(Error_prop, Error/100)]
    fixed_y_position <- max(plot_data$Proportion + plot_data$Error_prop, na.rm = TRUE) * 1.15
    fixed_y_position <- max(fixed_y_position, max(plot_data$Proportion) * 1.2)
    fixed_y_position <- min(fixed_y_position, 0.7)
    max_y <- 0.6
    plot_data[, `:=`(Sig_label, ifelse(P_value < 0.05, ifelse(P_value < 0.01, ifelse(P_value < 0.001, "***", "**"), "*"), 
        " ns"))]
    plot_data[, `:=`(Label, paste0("OR = ", round(OR, 2), Sig_label))]
    p <- ggplot(plot_data, aes(x = Region, y = Proportion, fill = Category)) + geom_col(position = position_dodge(0.8), width = 0.7, 
        color = "white", linewidth = 0.3) + geom_errorbar(aes(ymin = Proportion - Error_prop, ymax = Proportion + Error_prop), 
        position = position_dodge(0.8), width = 0.15, linewidth = 0.8, color = "#F1DD10", alpha = 0.8, na.rm = TRUE) + geom_errorbar(aes(ymin = Proportion - 
        Error_prop, ymax = Proportion - Error_prop), position = position_dodge(0.8), width = 0.3, linewidth = 0.8, color = "#F1DD10", 
        alpha = 0.8, na.rm = TRUE) + geom_errorbar(aes(ymin = Proportion + Error_prop, ymax = Proportion + Error_prop), position = position_dodge(0.8), 
        width = 0.3, linewidth = 0.8, color = "#F1DD10", alpha = 0.8, na.rm = TRUE) + geom_text(aes(x = as.numeric(Region) + 
        x_offset, y = fixed_y_position, label = Label), size = 3.5, angle = 45, hjust = 0.5, vjust = 0, color = ifelse(plot_data$Category == 
        "Allosteric", "#F4270C", ifelse(plot_data$Category == "Inhibit", "#75C2F6", "#F4AD0C"))) + theme_classic(base_size = 15) + 
        theme(plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 15)), axis.title = element_text(size = 14), 
            axis.text = element_text(size = 12, color = "black"), axis.text.x = element_text(size = 13), axis.line = element_line(color = "black", 
                linewidth = 0.5), axis.ticks = element_line(color = "black", linewidth = 0.5), legend.title = element_text(size = 13), 
            legend.text = element_text(size = 12), legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal", 
            panel.grid = element_blank(), plot.margin = margin(t = 40, r = 10, b = 10, l = 10)) + geom_hline(yintercept = 0, 
        color = "black", linewidth = 0.5) + labs(title = paste0("Enrichment of Allosteric Mutations - ", assay_name), y = "Proportion of all allosteric mutations", 
        fill = "Direction") + scale_y_continuous(limits = c(0, max_y), expand = expansion(mult = c(0, 0.02)), breaks = seq(0, 
        max_y, 0.1), labels = scales::percent) + scale_fill_manual(values = category_colors)
    ggsave(filename = file.path(output_dir, paste0("allosteric mutations Enrichment_plot_", assay_name, ".pdf")), plot = p, 
        width = 8, height = 6, dpi = 300)
    ggsave(filename = file.path(output_dir, paste0("allosteric mutations Enrichment_plot_", assay_name, ".png")), plot = p, 
        width = 8, height = 6, dpi = 300)
    return(p)
}

