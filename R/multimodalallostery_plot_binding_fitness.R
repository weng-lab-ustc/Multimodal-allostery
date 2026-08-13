#' Plot Binding Fitness.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input Input file path or input data object.
#' @param assay_sele Selected assay identifier.
#' @param anno Annotation data or annotation file path used by the analysis.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_binding_fitness <- function(input, assay_sele, anno) {
    input_abundance <- input[assay == "stab", ]
    input_abundance_single <- krasddpcams__nor_overlap_single_mt_fitness(input_abundance)
    input_binding <- input[assay == assay_sele, ]
    input_binding_single <- krasddpcams__nor_overlap_single_mt_fitness(input_binding)
    input_long <- rbind(input_abundance_single, input_binding_single)
    input_dc <- dcast(input_long, nt_seq + aa_seq + Nham_aa + AA_Pos1 + wtcodon1 ~ assay, value.var = c("nor_fitness_nooverlap", 
        "nor_fitness_nooverlap_sigma", "nor_gr_nooverlap", "nor_gr_nooverlap_sigma"), drop = TRUE)
    input_single_pos <- input_dc
    input_single_pos[, `:=`(position, AA_Pos1)]
    input_single_pos[, `:=`(WT_AA, wtcodon1)]
    anno_single <- merge(input_single_pos, anno, by.x = c("position", "WT_AA"), by.y = c("Pos_real", "codon"), all = TRUE)
    interface_col <- paste0("scHAmin_ligand_", assay_sele)
    anno_single[, `:=`(type_bs, "others")]
    anno_single[get(interface_col) < 5, `:=`(type_bs, "binding_interface")]
    p <- ggplot() + geom_point(data = anno_single[position > 1 & type_bs == "others", ], aes(x = nor_fitness_nooverlap_stab, 
        y = get(paste0("nor_fitness_nooverlap_", assay_sele))), color = "grey70", alpha = 0.4, size = 1.5) + geom_point(data = anno_single[position > 
        1 & type_bs == "binding_interface", ], aes(x = nor_fitness_nooverlap_stab, y = get(paste0("nor_fitness_nooverlap_", 
        assay_sele))), color = "#F4270C", alpha = 0.5, size = 2) + coord_cartesian(xlim = c(-1.3, 1.1), ylim = c(-1.3, 0.3)) + 
        theme_classic() + theme(text = element_text(size = 10), axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, 
        vjust = 0.5, hjust = 1), legend.position = "right", legend.key.height = unit(3.1, "mm"), legend.key.width = unit(3.1, 
        "mm"), plot.margin = margin(0, 0, 0, 0)) + labs(x = "Abundance Fitness", y = paste0(assay_sele, " Binding Fitness"))
    return(p)
}

