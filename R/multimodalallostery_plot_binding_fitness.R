#' plot binding fitness
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_binding_fitness()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_sele Input or option used by this workflow; see Usage for the exact interface.
#' @param anno Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_binding_fitness <- function (input, assay_sele, anno) 
{
    input_abundance <- input[assay == "stab", ]
    input_abundance_single <- krasddpcams::krasddpcams__nor_overlap_single_mt_fitness(input_abundance)
    input_binding <- input[assay == assay_sele, ]
    input_binding_single <- krasddpcams::krasddpcams__nor_overlap_single_mt_fitness(input_binding)
    input_long <- rbind(input_abundance_single, input_binding_single)
    input_dc <- data.table::dcast(input_long, nt_seq + aa_seq + Nham_aa + AA_Pos1 + wtcodon1 ~ assay, value.var = c("nor_fitness_nooverlap", 
        "nor_fitness_nooverlap_sigma", "nor_gr_nooverlap", "nor_gr_nooverlap_sigma"), drop = TRUE)
    input_single_pos <- input_dc
    input_single_pos[, `:=`(position, AA_Pos1)]
    input_single_pos[, `:=`(WT_AA, wtcodon1)]
    anno_single <- merge(input_single_pos, anno, by.x = c("position", "WT_AA"), by.y = c("Pos_real", "codon"), all = TRUE)
    interface_col <- paste0("scHAmin_ligand_", assay_sele)
    anno_single[, `:=`(type_bs, "others")]
    anno_single[get(interface_col) < 5, `:=`(type_bs, "binding_interface")]
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = anno_single[position > 1 & type_bs == "others", ], ggplot2::aes(x = nor_fitness_nooverlap_stab, 
        y = get(paste0("nor_fitness_nooverlap_", assay_sele))), color = "grey70", alpha = 0.40000000000000002, size = 1.5) + 
        ggplot2::geom_point(data = anno_single[position > 1 & type_bs == "binding_interface", ], ggplot2::aes(x = nor_fitness_nooverlap_stab, 
            y = get(paste0("nor_fitness_nooverlap_", assay_sele))), color = "#F4270C", alpha = 0.5, size = 2) + ggplot2::coord_cartesian(xlim = c(-1.3, 
        1.1000000000000001), ylim = c(-1.3, 0.29999999999999999)) + ggplot2::theme_classic() + ggplot2::theme(text = ggplot2::element_text(size = 10), 
        axis.text = ggplot2::element_text(size = 10), axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "right", legend.key.height = grid::unit(3.1000000000000001, "mm"), legend.key.width = grid::unit(3.1000000000000001, 
            "mm"), plot.margin = ggplot2::margin(0, 0, 0, 0)) + ggplot2::labs(x = "Abundance Fitness", y = paste0(assay_sele, 
        " Binding Fitness"), color = NULL)
    return(p)
}

