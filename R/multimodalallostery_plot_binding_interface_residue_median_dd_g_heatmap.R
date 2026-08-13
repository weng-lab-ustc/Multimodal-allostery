#' Plot Binding Interface Residue Median Dd g Heatmap.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param ddG_file Value supplied for `ddG_file`.
#' @param binding_sites Value supplied for `binding_sites`.
#' @param position_labels Value supplied for `position_labels`.
#' @param title Value supplied for `title`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_binding_interface_residue_median_dd_g_heatmap <- function(ddG_file, binding_sites, position_labels, title = "Binding Interface Residues - Median \u0394\u0394Gb") {
    ddG <- fread(ddG_file)
    ddG <- ddG[, `:=`(Pos_real, Pos + 1)]
    median_values <- ddG %>% filter(Pos_real %in% binding_sites) %>% group_by(Pos_real) %>% summarise(median_ddG = median(`mean_kcal/mol`, 
        na.rm = TRUE)) %>% ungroup() %>% mutate(residue_label = factor(position_labels[as.character(Pos_real)], levels = position_labels), 
        ) %>% arrange(residue_label)
    p <- ggplot(median_values, aes(x = 1, y = residue_label, fill = median_ddG)) + geom_tile(color = "white", linewidth = 0.1) + 
        scale_fill_gradient2(low = "#1B38A6", mid = "gray", high = "#F4270C", midpoint = 0, limits = c(-1, 2.5), name = expression(Delta * 
            Delta * "Gb (kcal/mol)")) + labs(title = title, y = "Residue") + theme_minimal() + theme(text = element_text(size = 8), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, 
            size = 8), axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 8), legend.position = "right", 
        panel.border = element_rect(color = "gray90", fill = NA, linewidth = 0.8), legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8)) + scale_y_discrete(limits = rev(levels(median_values$residue_label))) + coord_fixed(ratio = 0.4)
    return(p)
}

