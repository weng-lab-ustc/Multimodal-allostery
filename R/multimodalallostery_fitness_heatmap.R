#' Fitness Heatmap.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input Input file path or input data object.
#' @param wt_aa Value supplied for `wt_aa`.
#' @param title Value supplied for `title`.
#' @param legend_limits Value supplied for `legend_limits`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_fitness_heatmap <- function(input, wt_aa, title = "fitness", legend_limits = c(-1.5, 1)) {
    aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
    num <- nchar(wt_aa) + 1
    input_single <- input
    input_single[, `:=`(position, AA_Pos1)]
    input_single[, `:=`(WT_AA, wtcodon1)]
    heatmap_tool_fitness <- data.table(wtcodon1 = rep(unlist(strsplit(wt_aa, "")), each = 21), position = rep(2:num, each = 21), 
        codon1 = c(unlist(aa_list), "*"))
    heatmap_tool_fitness_anno_single <- merge(input_single, heatmap_tool_fitness, by = c("wtcodon1", "position", "codon1"), 
        all = TRUE)
    heatmap_tool_fitness_anno_single <- within(heatmap_tool_fitness_anno_single, codon1 <- factor(codon1, levels = c("*", 
        "D", "E", "R", "H", "K", "S", "T", "N", "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "W", "Y")))
    heatmap_tool_fitness_anno_single[wtcodon1 == codon1, `:=`(nor_fitness_nooverlap, 0)]
    ggplot() + theme_classic() + geom_tile(data = heatmap_tool_fitness_anno_single[position > 1, ], aes(x = position, y = codon1, 
        fill = nor_fitness_nooverlap)) + scale_x_discrete(limits = c(2:num), labels = c(2:num)) + theme(axis.text.x = element_text(size = 8, 
        vjust = 0.5, hjust = 0.5, color = c(NA, NA, NA, rep(c("black", NA, NA, NA, NA), nchar(wt_aa)%/%5)))) + scale_fill_gradient2(limits = legend_limits, 
        low = "#F4270C", mid = "gray", high = "#1B38A6", name = "Fitness", midpoint = 0, na.value = "white", guide = guide_colorbar(title.position = "top", 
            title.hjust = 0.5)) + ylab("Mutant aa") + ggtitle(title) + labs() + geom_text(data = heatmap_tool_fitness_anno_single[position > 
        1 & wtcodon1 == codon1, ], aes(x = position, y = codon1), label = "-", size = 3) + theme(text = element_text(size = 8), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), legend.position = c(1, 1.38), title = element_text(size = 8), 
        legend.justification = c(1, 1), legend.direction = "horizontal", legend.text = element_text(size = 8), axis.title.x = element_text(size = 8, 
            face = "plain"), axis.title.y = element_text(size = 8, face = "plain"), axis.text.x = element_text(size = 8, 
            angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(family = "Courier", angle = 90, size = 9.5, vjust = 0.5, 
            hjust = 0.5, margin = margin(0, -0.5, 0, 0, "mm")), legend.key.height = unit(3.1, "mm"), legend.key.width = unit(4, 
            "mm"), legend.key.size = unit(1, "mm"), plot.margin = margin(0, -0, 0, 0)) + coord_fixed()
}

