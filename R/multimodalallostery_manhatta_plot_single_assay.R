#' Manhatta Plot Single Assay.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input Input file path or input data object.
#' @param assay_sele Selected assay identifier.
#' @param anno Annotation data or annotation file path used by the analysis.
#' @param rects_sheet Value supplied for `rects_sheet`.
#' @param rects_alpha Value supplied for `rects_alpha`.
#' @param wt_aa Value supplied for `wt_aa`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_manhatta_plot_single_assay <- function(input, assay_sele, anno, rects_sheet, rects_alpha, wt_aa) {
    ddG <- fread(input)
    ddG[, `:=`(Pos_real, Pos_ref + 1)]
    ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
    ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
    ddG[, `:=`(mt, paste0(wt_codon, Pos_real, mt_codon))]
    aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
    heatmap_tool <- data.table(wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20), Pos_real = rep(2:188, each = 20), mt_codon = rep(aa_list, 
        times = length(strsplit(wt_aa, "")[[1]])))
    ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
    ddG[, `:=`(Pos, Pos_real)]
    output <- ddG[Pos_real > 1, .(mean = sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = TRUE)/sum(1/.SD[[2]]^2, na.rm = TRUE)), .SDcols = c("mean_kcal/mol", 
        "std_kcal/mol"), by = "Pos_real"]
    output_sigma <- ddG[Pos_real > 1, .(sigma = sqrt(1/sum(1/.SD[[2]]^2, na.rm = TRUE))), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), 
        by = "Pos_real"]
    weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
    weighted_mean_ddG[, `:=`(Pos, Pos_real)]
    anno <- fread(anno)
    data_plot <- merge(weighted_mean_ddG, anno, by = "Pos", all = TRUE)
    data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5, `:=`(binding_type, "binding site")]
    data_plot[, `:=`(binding_type_gtp_included, binding_type)]
    data_plot[get(paste0("GXPMG_scHAmin_ligand_", assay_sele)) < 5, `:=`(binding_type_gtp_included, "GTP binding site")]
    reg_threshold <- data_plot[binding_type == "binding site", sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = TRUE)/sum(1/.SD[[2]]^2, 
        na.rm = TRUE), .SDcols = c("mean", "sigma")]
    data_plot[, `:=`(site_type, "Reminder")]
    data_plot[binding_type_gtp_included == "binding site", `:=`(site_type, "Binding interface site")]
    data_plot[binding_type_gtp_included == "GTP binding site", `:=`(site_type, "GTP binding interface site")]
    data_plot_mutation1 <- merge(ddG, data_plot[, .(Pos, site_type)], by = "Pos", all.x = TRUE)
    data_plot_mutation <- data_plot_mutation1[Pos > 1 & !is.na(id)]
    data_plot_mutation[, `:=`(mutation_type, "Reminder")]
    data_plot_mutation[, `:=`(allosteric_mutation, p.adjust(krasddpcams__pvalue(abs(mean) - reg_threshold, std), method = "BH") < 
        0.05 & (abs(mean) - reg_threshold) > 0)]
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
    p <- ggplot() + geom_rect(data = rects_sheet, aes(ymin = -3.5, ymax = 3.5, xmin = xstart - 0.5, xmax = xend + 0.5), fill = "#75C2F6", 
        alpha = 0.06) + geom_rect(data = rects_alpha, aes(ymin = -3.5, ymax = 3.5, xmin = xstart - 0.5, xmax = xend + 0.5), 
        fill = "#C68EFD", alpha = 0.06) + geom_point(data = data_plot_mutation, aes(x = Pos_real, y = `mean_kcal/mol`, color = mutation_type), 
        size = 1) + scale_color_manual(values = c(alpha("#F4270C", 1), alpha("#FFB0A5", 1), alpha("#1B38A6", 1), alpha("#75C2F6", 
        0.6), alpha("#F4AD0C", 1), alpha("gray", 0.8))) + geom_hline(yintercept = 0, linetype = 2) + scale_x_continuous(expand = c(1/188, 
        11/188)) + ylab(paste0("Binding \u0394\u0394G (", assay_sele, ") (kcal/mol)")) + xlab("Amino acid position") + labs() + annotate("text", 
        x = (rects_sheet$xstart[1] + rects_sheet$xend[1])/2, y = 3.1, label = "strand", size = 3.5, vjust = 0.5, hjust = 1, 
        angle = 90) + annotate("text", x = (rects_alpha$xstart[1] + rects_alpha$xend[1])/2, y = 3.1, label = "helix", size = 3.5, 
        vjust = 0.5, hjust = 1, angle = 90) + theme_classic2() + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
        text = element_text(size = 12), legend.position = "none", legend.text = element_text(size = 12), strip.background = element_rect(colour = "black", 
            fill = "white")) + coord_fixed(ratio = 10, xlim = c(-0.5, 190), ylim = c(-1.5, 3.5))
    return(list(plot = p, mutation_count = data_plot_mutation[, .N, by = mutation_type]))
}

