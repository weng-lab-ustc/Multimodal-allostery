#' Plot Dd g Correlation.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data1 Value supplied for `data1`.
#' @param data2 Value supplied for `data2`.
#' @param data1_name Value supplied for `data1_name`.
#' @param data2_name Value supplied for `data2_name`.
#' @param output_file Optional output file path.
#' @param limits Value supplied for `limits`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_dd_g_correlation <- function(data1, data2, data1_name = "Dataset 1", data2_name = "Dataset 2", output_file = NULL, limits = c(-1.6, 
    2.8)) {
    points_per_mm <- 72.27/25.4
    data1_dt <- as.data.table(data1)
    data2_dt <- as.data.table(data2)
    cat("Checking for duplicate positions...\n")
    cat("Data1 duplicate positions:", sum(duplicated(data1_dt$Pos)), "\n")
    cat("Data2 duplicate positions:", sum(duplicated(data2_dt$Pos)), "\n")
    merged_data <- merge(data1_dt, data2_dt, by = c("Pos", "id"), suffixes = c("_1", "_2"))
    x_values <- merged_data[["mean_kcal/mol_1"]]
    y_values <- merged_data[["mean_kcal/mol_2"]]
    cor_test <- cor.test(x_values, y_values, method = "pearson")
    r_value <- round(cor_test$estimate, 3)
    p_value <- cor_test$p.value
    p_text <- ifelse(p_value < 0.001, "p < 0.001", ifelse(p_value < 0.01, "p < 0.01", ifelse(p_value < 0.05, "p < 0.05", 
        paste0("p = ", round(p_value, 3)))))
    p <- ggplot(merged_data, aes(x = x_values, y = y_values)) + geom_point(alpha = 0.35, size = 0.6, color = "#75C2F6") + 
        annotate("text", x = -Inf, y = Inf, label = paste0("R = ", r_value, "\n", p_text), hjust = -0.1, vjust = 1.5, size = 8/points_per_mm, 
            color = "black") + scale_x_continuous(limits = limits, breaks = scales::pretty_breaks(n = 6)) + scale_y_continuous(limits = limits, 
        breaks = scales::pretty_breaks(n = 6)) + labs(x = paste(data1_name, "\u0394\u0394G (kcal/mol)"), y = paste(data2_name, "\u0394\u0394G (kcal/mol)"), 
        title = "\u0394\u0394G Correlation Analysis") + theme_classic() + theme(text = element_text(size = 8), axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8, color = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5), plot.title = element_text(size = 8, hjust = 0.5), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4)) + 
        coord_fixed(ratio = 1)
    cat("Correlation Analysis:\n")
    cat("Number of positions:", nrow(merged_data), "\n")
    cat("Pearson R:", r_value, "\n")
    cat("P-value:", p_value, "\n")
    cat("P-value display:", p_text, "\n")
    if (!is.null(output_file)) {
        ggsave(output_file, p, device = cairo_pdf, width = 4, height = 4, units = "in", dpi = 300)
        cat("Plot saved to:", output_file, "\n")
    }
    return(p)
}

