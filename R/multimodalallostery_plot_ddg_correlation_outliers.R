#' Plot Ddg Correlation Outliers.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param analysis Value supplied for `analysis`.
#' @param x_label Value supplied for `x_label`.
#' @param y_label Value supplied for `y_label`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_ddg_correlation_outliers <- function(analysis, x_label, y_label) {
    data <- analysis$data
    outliers <- analysis$outliers
    x_var <- analysis$x_var
    y_var <- analysis$y_var
    r_value <- analysis$correlation$estimate
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])) + ggplot2::geom_smooth(method = "lm", 
        se = FALSE, color = "black", linewidth = 0.8, alpha = 0.5) + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[[y_var]] - 
        `std_kcal/mol.y`, ymax = .data[[y_var]] + `std_kcal/mol.y`), width = 0.05, color = "grey70", alpha = 0.3) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[x_var]] - 
        `std_kcal/mol.x`, xmax = .data[[x_var]] + `std_kcal/mol.x`), height = 0.05, color = "grey70", alpha = 0.3) + ggplot2::geom_point(color = "#75C2F6", 
        alpha = 0.8, size = 2) + ggplot2::geom_errorbar(data = outliers, ggplot2::aes(ymin = .data[[y_var]] - `std_kcal/mol.y`, 
        ymax = .data[[y_var]] + `std_kcal/mol.y`), width = 0.08, color = "#FF6A56", alpha = 0.6, linewidth = 0.8) + ggplot2::geom_errorbarh(data = outliers, 
        ggplot2::aes(xmin = .data[[x_var]] - `std_kcal/mol.x`, xmax = .data[[x_var]] + `std_kcal/mol.x`), height = 0.08, 
        color = "#FF6A56", alpha = 0.6, linewidth = 0.8) + ggplot2::geom_point(data = outliers, color = "#FF6A56", size = 2.5) + 
        ggrepel::geom_text_repel(data = outliers, ggplot2::aes(label = paste0(mt, " (\u8796=", round(abs_residual, 2), ")")), 
            color = "#FF6A56", size = 2.5, max.overlaps = Inf, box.padding = 0.7, point.padding = 0.3, segment.color = "#FF6A56", 
            segment.alpha = 0.3) + ggplot2::annotate("text", x = min(data[[x_var]], na.rm = TRUE), y = max(data[[y_var]], 
        na.rm = TRUE), label = paste("Pearson r =", round(r_value, 3), "\nTop", analysis$num_outliers, "deviating points labeled\nError bars show \u53641 SD"), 
        hjust = 0, vjust = 1, size = 2.5) + ggplot2::labs(title = "Correlation Analysis with Error Bars - Most Deviating Points", 
        x = x_label, y = y_label, caption = "Error bars represent standard deviation.\nPoints with large residuals that cannot be explained by measurement error are likely true biological outliers.") + 
        ggplot2::theme_minimal() + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.line = ggplot2::element_line(linewidth = 1), 
        axis.ticks = ggplot2::element_line(linewidth = 1), axis.text = ggplot2::element_text(size = 8), axis.text.x = ggplot2::element_text(angle = 90, 
            vjust = 0.5, hjust = 1), axis.title = ggplot2::element_text(size = 8), plot.title = ggplot2::element_text(size = 8), 
        plot.caption = ggplot2::element_text(size = 8, color = "grey50", hjust = 0), legend.text = ggplot2::element_text(size = 8), 
        legend.title = ggplot2::element_text(size = 8))
    return(p)
}

