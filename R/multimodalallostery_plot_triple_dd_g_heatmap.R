#' Plot Triple Dd g Heatmap.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param ddG_K13 Value supplied for `ddG_K13`.
#' @param ddG_K19 Value supplied for `ddG_K19`.
#' @param ddG_RAF1 Value supplied for `ddG_RAF1`.
#' @param anno Annotation data or annotation file path used by the analysis.
#' @param wt_aa Value supplied for `wt_aa`.
#' @param colour_scheme Value supplied for `colour_scheme`.
#' @param allosteric_sites_list Value supplied for `allosteric_sites_list`.
#' @param legend_limits Value supplied for `legend_limits`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_triple_dd_g_heatmap <- function(ddG_K13, ddG_K19, ddG_RAF1, anno, wt_aa, colour_scheme, allosteric_sites_list = NULL, 
    legend_limits = c(-1.3, 3)) {
    aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
    heatmap_tool <- data.table(wt_codon = rep(unlist(strsplit(wt_aa, "")), each = 20), Pos_real = rep(2:188, each = 20), 
        mt_codon = unlist(aa_list))
    ddG1 <- krasddpcams__read_ddG(ddG_K13, "K13")
    input1_heatmap <- merge(ddG1, heatmap_tool, by = c("wt_codon", "Pos_real", "mt_codon"), all = T)
    ddG2 <- krasddpcams__read_ddG(ddG_K19, "K19")
    input2_heatmap <- merge(ddG2, heatmap_tool, by = c("wt_codon", "Pos_real", "mt_codon"), all = T)
    ddG3 <- krasddpcams__read_ddG(ddG_RAF1, "RAF1")
    input3_heatmap <- merge(ddG3, heatmap_tool, by = c("wt_codon", "Pos_real", "mt_codon"), all = T)
    ddG <- rbind(input1_heatmap, input2_heatmap, input3_heatmap)
    ddG[wt_codon == mt_codon, `:=`(`mean_kcal/mol`, 0)]
    output <- merge(ddG, anno, by.x = "Pos_real", by.y = "Pos", all = T)
    if (!is.null(allosteric_sites_list)) {
        all_sites <- unique(c(allosteric_sites_list$K13, allosteric_sites_list$K19, allosteric_sites_list$RAF1))
        output <- output[Pos_real %in% all_sites, ]
    }
    output[, `:=`(Pos_real, factor(Pos_real, levels = sort(unique(Pos_real))))]
    output <- within(output, mt_codon <- factor(mt_codon, levels = c("D", "E", "R", "H", "K", "S", "T", "N", "Q", "C", "G", 
        "P", "A", "V", "I", "L", "M", "F", "W", "Y")))
    output <- within(output, assay <- factor(assay, levels = c("K13", "K19", "RAF1")))
    output[, `:=`(wtcodon_pos, paste0(wt_codon, Pos_real))]
    output[, `:=`(Pos_num, as.numeric(as.character(Pos_real)))]
    output <- output[order(Pos_num), ]
    output[, `:=`(wtcodon_pos, factor(wtcodon_pos, levels = unique(wtcodon_pos)))]
    p <- ggplot2::ggplot() + geom_tile(data = output, aes(x = assay, y = mt_codon, fill = `mean_kcal/mol`)) + geom_text(data = output[Pos_num > 
        1 & wt_codon == mt_codon, ], aes(x = assay, y = mt_codon, label = "-"), size = 6 * 5/14) + scale_fill_gradient2(limits = legend_limits, 
        low = colour_scheme[["blue"]], mid = "gray", high = colour_scheme[["red"]], na.value = "white") + facet_wrap(~wtcodon_pos, 
        nrow = 2) + ylab("Mutant AA") + xlab("Binding partners") + labs(fill = expression(Delta * Delta * G ~ "(kcal/mol)")) + 
        theme_bw() + theme(text = element_text(size = 15, family = "Arial"), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "bottom", strip.background = element_rect(colour = "white", fill = "white"), axis.text.x = element_text(angle = 90, 
            vjust = 0.5, hjust = 1), axis.text.y = element_text(margin = margin(0, 5, 0, 0)), panel.spacing.y = unit(3, "mm"), 
        panel.spacing.x = unit(1, "mm")) + coord_fixed()
    return(p)
}

