#' manhatta plot single assay
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_manhatta_plot_single_assay()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_sele Input or option used by this workflow; see Usage for the exact interface.
#' @param anno Input or option used by this workflow; see Usage for the exact interface.
#' @param rects_sheet Input or option used by this workflow; see Usage for the exact interface.
#' @param rects_alpha Input or option used by this workflow; see Usage for the exact interface.
#' @param wt_aa Input or option used by this workflow; see Usage for the exact interface.
#' @param threshold Input or option used by this workflow; see Usage for the exact interface.
#' @param print_summary Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_manhatta_plot_single_assay <- function (input, assay_sele, anno, rects_sheet, rects_alpha, wt_aa, threshold = 0.40000000000000002, print_summary = TRUE) 
{
    ddG <- .ma_read_table(input)
    ddG[, `:=`(Pos_real, Pos_ref + 1)]
    ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
    ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
    ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
    aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
    heatmap_tool <- data.table::data.table(wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20), Pos_real = rep(2:188, each = 20), 
        mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]])))
    ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
    ddG[, `:=`(Pos, Pos_real)]
    output <- ddG[Pos_real > 1, .(mean = sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = TRUE)/sum(1/.SD[[2]]^2, na.rm = TRUE)), .SDcols = c("mean_kcal/mol", 
        "std_kcal/mol"), by = "Pos_real"]
    output_sigma <- ddG[Pos_real > 1, .(sigma = sqrt(1/sum(1/.SD[[2]]^2, na.rm = TRUE))), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), 
        by = "Pos_real"]
    weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
    weighted_mean_ddG[, `:=`(Pos, Pos_real)]
    anno <- .ma_read_table(anno)
    data_plot <- merge(weighted_mean_ddG, anno, by = "Pos", all = TRUE)
    data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, `:=`(binding_type, "binding site")]
    data_plot[, `:=`(binding_type_gtp_included, binding_type)]
    data_plot[get(paste0("GXPMG_scHAmin_ligand_", assay_sele)) < 5, `:=`(binding_type_gtp_included, "GTP binding site")]
    reg_threshold <- threshold
    data_plot[, `:=`(site_type, "Reminder")]
    data_plot[binding_type_gtp_included == "binding site", `:=`(site_type, "Binding interface site")]
    data_plot[binding_type_gtp_included == "GTP binding site", `:=`(site_type, "GTP binding interface site")]
    data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
    data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
    data_plot_mutation[, `:=`(mutation_type, "Reminder")]
    data_plot_mutation[, `:=`(allosteric_mutation, stats::p.adjust(krasddpcams::krasddpcams__pvalue(abs(`mean_kcal/mol`) - 
        reg_threshold, `std_kcal/mol`), method = "BH") < 0.050000000000000003 & (abs(`mean_kcal/mol`) - reg_threshold) > 
        0)]
    data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] & allosteric_mutation == TRUE, `:=`(mutation_type, 
        "Orthosteric site huge differences")]
    data_plot_mutation[Pos %in% data_plot[site_type == "Binding interface site", Pos] & allosteric_mutation == FALSE, `:=`(mutation_type, 
        "Orthosteric site small differences")]
    data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] & allosteric_mutation == TRUE, 
        `:=`(mutation_type, "GTP binding allosteric mutation")]
    data_plot_mutation[Pos %in% data_plot[site_type == "GTP binding interface site", Pos] & allosteric_mutation == FALSE, 
        `:=`(mutation_type, "GTP binding other mutation")]
    data_plot_mutation[!site_type %in% c("GTP binding interface site", "Binding interface site") & allosteric_mutation == 
        TRUE, `:=`(mutation_type, "Allosteric mutation")]
    data_plot_mutation[!site_type %in% c("GTP binding interface site", "Binding interface site") & allosteric_mutation == 
        FALSE, `:=`(mutation_type, "Other mutation")]
    data_plot_mutation <- within(data_plot_mutation, mutation_type <- factor(mutation_type, levels = c("Orthosteric site huge differences", 
        "Orthosteric site small differences", "GTP binding allosteric mutation", "GTP binding other mutation", "Allosteric mutation", 
        "Other mutation")))
    if (print_summary) {
        cat("\n", paste0(rep("=", 80), collapse = ""), "\n")
        cat("MANHATTAN PLOT SUMMARY FOR:", assay_sele, "\n")
        cat(paste0(rep("=", 80), collapse = ""), "\n")
        cat("Fixed threshold used:", threshold, "kcal/mol\n")
        cat("Total mutations analyzed:", nrow(data_plot_mutation), "\n\n")
        mutation_counts <- data_plot_mutation[, .N, by = mutation_type]
        data.table::setorder(mutation_counts, -N)
        cat("--- Mutation Type Distribution ---\n")
        for (i in 1:nrow(mutation_counts)) {
            cat(sprintf("  %-35s: %3d\n", mutation_counts$mutation_type[i], mutation_counts$N[i]))
        }
        cat(paste0(rep("=", 80), collapse = ""), "\n\n")
    }
    p <- ggplot2::ggplot() + ggplot2::geom_rect(data = rects_sheet, ggplot2::aes(ymin = -3.5, ymax = 3.5, xmin = xstart - 
        0.5, xmax = xend + 0.5), fill = "#75C2F6", alpha = 0.059999999999999998) + ggplot2::geom_rect(data = rects_alpha, 
        ggplot2::aes(ymin = -3.5, ymax = 3.5, xmin = xstart - 0.5, xmax = xend + 0.5), fill = "#C68EFD", alpha = 0.059999999999999998) + 
        ggplot2::geom_point(data = data_plot_mutation, ggplot2::aes(x = Pos_real, y = `mean_kcal/mol`, color = mutation_type), 
            size = 1) + ggplot2::scale_color_manual(values = c(scales::alpha("#F4270C", 1), scales::alpha("#FFB0A5", 1), 
        scales::alpha("#1B38A6", 1), scales::alpha("#75C2F6", 0.59999999999999998), scales::alpha("#F4AD0C", 1), scales::alpha("gray", 
            0.80000000000000004))) + ggplot2::geom_hline(yintercept = 0, linetype = 2) + ggplot2::geom_hline(yintercept = threshold, 
        linetype = 3, color = "red", alpha = 0.5) + ggplot2::geom_hline(yintercept = -threshold, linetype = 3, color = "red", 
        alpha = 0.5) + ggplot2::scale_x_continuous(expand = c(1/188, 11/188)) + ggplot2::ylab(paste0("Binding \u0394\u0394G (", assay_sele, 
        ") (kcal/mol)")) + ggplot2::xlab("Amino acid position") + ggplot2::labs(color = NULL) + ggplot2::annotate("text", 
        x = (rects_sheet$xstart[1] + rects_sheet$xend[1])/2, y = 3.1000000000000001, label = "strand", size = 3.5, vjust = 0.5, 
        hjust = 1, angle = 90) + ggplot2::annotate("text", x = (rects_alpha$xstart[1] + rects_alpha$xend[1])/2, y = 3.1000000000000001, 
        label = "helix", size = 3.5, vjust = 0.5, hjust = 1, angle = 90) + ggpubr::theme_classic2() + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12), 
        axis.text.y = ggplot2::element_text(size = 12), text = ggplot2::element_text(size = 12), legend.position = "none", 
        legend.text = ggplot2::element_text(size = 12), strip.background = ggplot2::element_rect(colour = "black", fill = "white")) + 
        ggplot2::coord_fixed(ratio = 10, xlim = c(-0.5, 190), ylim = c(-1.5, 3.5))
    return(list(plot = p, mutation_count = data_plot_mutation[, .N, by = mutation_type], threshold_used = threshold, allosteric_count = data_plot_mutation[allosteric_mutation == 
        TRUE, .N], total_mutations = nrow(data_plot_mutation), data = data_plot_mutation))
}

