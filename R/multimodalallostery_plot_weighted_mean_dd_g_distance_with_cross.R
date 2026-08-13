#' Plot Weighted Mean Dd g Distance with Cross.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param ddG_files Value supplied for `ddG_files`.
#' @param assays Value supplied for `assays`.
#' @param anno_file Value supplied for `anno_file`.
#' @param x_intercept Value supplied for `x_intercept`.
#' @param output_file Optional output file path.
#' @param base_font_size Value supplied for `base_font_size`.
#' @param point_size Value supplied for `point_size`.
#' @param text_repel_size Value supplied for `text_repel_size`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_weighted_mean_dd_g_distance_with_cross <- function(ddG_files, assays, anno_file, x_intercept = 5, output_file = NULL, 
    base_font_size = 12, point_size = 0.5, text_repel_size = 3) {
    rect_input <- data.frame(xstart = c(3, 15, 38, 51, 67, 77, 87, 109, 127, 139, 148), xend = c(9, 24, 44, 57, 73, 84, 104, 
        115, 136, 143, 166), col = c("b1", "a1", "b2", "b3", "a2", "b4", "a3", "b5", "a4", "b6", "a5"))
    anno <- fread(anno_file)
    anno[, `:=`(Pos_real, Pos)]
    weighted_mean_ddG <- list()
    for (i in seq_along(assays)) {
        weighted_mean_ddG[[assays[i]]] <- krasddpcams__get_weighted_mean_abs_ddG_mutcount(ddG = ddG_files[i], assay_sele = assays[i])
    }
    data_plot <- data.table()
    for (assayi in assays) {
        data_plot_assayi <- merge(weighted_mean_ddG[[assayi]], anno, by = "Pos_real", all = TRUE)
        data_plot_assayi[, `:=`(binding_type, "allosteric site")]
        data_plot_assayi[get(paste0("scHAmin_ligand_", assayi)) < x_intercept, `:=`(binding_type, "binding site")]
        data_plot_assayi[, `:=`(binding_type_gtp_included, binding_type)]
        data_plot_assayi[GXPMG_scHAmin_ligand_RAF1 < x_intercept, `:=`(binding_type_gtp_included, "GTP binding site")]
        data_plot <- rbind(data_plot, data_plot_assayi)
    }
    reg_threshold <- data_plot[binding_type == "binding site", sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = TRUE)/sum(1/.SD[[2]]^2, 
        na.rm = TRUE), .SDcols = c("mean", "sigma")]
    print(paste("Regression threshold (based on all binding sites):", reg_threshold))
    data_plot[, `:=`(site_type, "Reminder")]
    data_plot[binding_type_gtp_included == "binding site", `:=`(site_type, "Binding interface site")]
    data_plot[binding_type_gtp_included == "GTP binding site", `:=`(site_type, "Other GTP pocket site")]
    data_plot[binding_type_gtp_included == "GTP binding site" & mean > reg_threshold & binding_type != "binding site" & count > 
        9.5, `:=`(site_type, "Allosteric GTP pocket site")]
    data_plot[binding_type_gtp_included == "allosteric site" & mean > reg_threshold & count > 9.5, `:=`(site_type, "Major allosteric site")]
    data_plot[, `:=`(colors_type, "others")]
    rects_dt <- as.data.table(rect_input)
    for (b in c("b1", "b2", "b3", "b4", "b5", "b6")) {
        data_plot[Pos_real >= rects_dt[col == b, xstart] & Pos_real <= rects_dt[col == b, xend], `:=`(colors_type, b)]
    }
    for (a in c("a1", "a2", "a3", "a4", "a5")) {
        data_plot[Pos_real >= rects_dt[col == a, xstart] & Pos_real <= rects_dt[col == a, xend], `:=`(colors_type, a)]
    }
    data_plot[, `:=`(shape, "others")]
    data_plot[colors_type %chin% c("b1", "b2", "b3", "b4", "b5", "b6"), `:=`(shape, "beta strand")]
    data_plot[colors_type %chin% c("a1", "a2", "a3", "a4", "a5"), `:=`(shape, "alpha helix")]
    data_plot <- data_plot[Pos_real > 1 & count > 9.5, ]
    data_plot <- within(data_plot, site_type <- factor(site_type, levels = c("Binding interface site", "Allosteric GTP pocket site", 
        "Other GTP pocket site", "Major allosteric site", "Reminder")))
    data_plot <- within(data_plot, assay <- factor(assay, levels = assays))
    allosteric_list <- list()
    for (assayi in assays) {
        data_plot[assay == assayi, `:=`(distance_bp, get(paste0("scHAmin_ligand_", assayi)))]
        allosteric_list[[assayi]] <- list()
        allosteric_list[[assayi]][["Binding interface site"]] <- data_plot[binding_type == "binding site" & assay == assayi, 
            Pos_real]
        allosteric_list[[assayi]][["Allosteric GTP pocket site"]] <- data_plot[site_type == "Allosteric GTP pocket site" & 
            assay == assayi, Pos_real]
        allosteric_list[[assayi]][["Other GTP pocket site"]] <- data_plot[site_type == "Other GTP pocket site" & assay == 
            assayi, Pos_real]
        allosteric_list[[assayi]][["Major allosteric site"]] <- data_plot[site_type == "Major allosteric site" & assay == 
            assayi, Pos_real]
    }
    protein_names <- assays
    for (current_protein in protein_names) {
        self_allo <- unique(c(allosteric_list[[current_protein]][["Allosteric GTP pocket site"]], allosteric_list[[current_protein]][["Major allosteric site"]]))
        other_proteins <- setdiff(protein_names, current_protein)
        other_hotspots <- data.table()
        for (other in other_proteins) {
            other_allo <- unique(c(allosteric_list[[other]][["Allosteric GTP pocket site"]], allosteric_list[[other]][["Major allosteric site"]]))
            if (length(other_allo) > 0) {
                other_hotspots <- rbind(other_hotspots, data.table(Pos_real = other_allo, source_protein = other))
            }
        }
        other_hotspots <- unique(other_hotspots)
        other_hotspots <- other_hotspots[!(Pos_real %in% self_allo), ]
        data_plot[assay == current_protein, `:=`(cross_hotspot_source, NA_character_)]
        data_plot[assay == current_protein, `:=`(label_color, NA_character_)]
        for (i in 1:nrow(other_hotspots)) {
            pos <- other_hotspots$Pos_real[i]
            source <- other_hotspots$source_protein[i]
            data_plot[assay == current_protein & Pos_real == pos, `:=`(cross_hotspot_source, source)]
            if (source == "K13") {
                data_plot[assay == current_protein & Pos_real == pos, `:=`(label_color, "#F1DD10")]
            }
            else if (source == "K19") {
                data_plot[assay == current_protein & Pos_real == pos, `:=`(label_color, "#C68EFD")]
            }
            else if (source == "RAF1") {
                data_plot[assay == current_protein & Pos_real == pos, `:=`(label_color, "#09B636")]
            }
        }
    }
    data_plot[, `:=`(point_color, as.character(site_type))]
    data_plot[!is.na(cross_hotspot_source), `:=`(point_color, paste0("Cross_", cross_hotspot_source))]
    color_values <- c(`Binding interface site` = "#F4270C", `Allosteric GTP pocket site` = "#1B38A6", `Major allosteric site` = "#F4AD0C", 
        `Other GTP pocket site` = "#75C2F6", Reminder = "gray", Cross_K13 = "#F1DD10", Cross_K19 = "#C68EFD", Cross_RAF1 = "#09B636")
    color_labels <- c(`Binding interface site` = "Binding interface", `Allosteric GTP pocket site` = "Allosteric GTP pocket (self)", 
        `Major allosteric site` = "Major allosteric (self)", `Other GTP pocket site` = "Other GTP pocket", Reminder = "Others", 
        Cross_K13 = "K13 allosteric hotspot", Cross_K19 = "K19 allosteric hotspot", Cross_RAF1 = "RAF1 allosteric hotspot")
    multimodalallostery_print_allosteric_statistics(data_plot, allosteric_list, assays, reg_threshold)
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("\u3010Cross-referenced hotspot annotation information (derived from other proteins and not constituting an allosteric hotspot of the protein itself)\u3011\n")
    cat(rep("-", 40), "\n", sep = "")
    for (current_protein in protein_names) {
        cross_data <- data_plot[assay == current_protein & !is.na(cross_hotspot_source), ]
        if (nrow(cross_data) > 0) {
            cat(sprintf("\nOther protein hotspots (non-autosteric hotspots) marked in %s:\n", current_protein))
            for (i in 1:nrow(cross_data)) {
                cat(sprintf("  Residue %d (from %s): mean ddG = %.3f, distance = %.2f \u00C5, point color = %s\n", cross_data$Pos_real[i], 
                  cross_data$cross_hotspot_source[i], cross_data$mean[i], cross_data$distance_bp[i], cross_data$cross_hotspot_source[i]))
            }
        }
        else {
            cat(sprintf("\nThere are no other protein hotspots requiring additional annotation in %s.\n", current_protein))
        }
    }
    cat(rep("=", 80), "\n\n", sep = "")
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = data_plot, mapping = aes(x = distance_bp, y = mean, color = point_color, 
        shape = as.factor(shape)), size = point_size) + ggplot2::geom_pointrange(data = data_plot, aes(x = distance_bp, y = mean, 
        color = point_color, ymin = mean - sigma, ymax = mean + sigma, shape = as.factor(shape)), size = point_size) + ggplot2::geom_hline(yintercept = reg_threshold, 
        linetype = 2, linewidth = 0.3) + ggplot2::geom_vline(xintercept = x_intercept, linetype = 2, linewidth = 0.3) + ggplot2::geom_hline(yintercept = 0, 
        linetype = "solid", linewidth = 0.3) + ggplot2::geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.3) + 
        ggrepel::geom_text_repel(data = data_plot[site_type == "Major allosteric site", ], aes(x = distance_bp, y = mean, 
            label = Pos_real), nudge_y = 0.05, color = "#F4AD0C", size = text_repel_size, fontface = "bold") + ggrepel::geom_text_repel(data = data_plot[site_type == 
        "Allosteric GTP pocket site", ], aes(x = distance_bp, y = mean, label = Pos_real), nudge_y = 0.05, color = "#1B38A6", 
        size = text_repel_size, fontface = "bold") + ggrepel::geom_text_repel(data = data_plot[!is.na(cross_hotspot_source) & 
        cross_hotspot_source == "K13", ], aes(x = distance_bp, y = mean, label = paste0(Pos_real, " (K13)")), nudge_y = 0.08, 
        size = text_repel_size * 0.8, fontface = "bold", color = "#F1DD10") + ggrepel::geom_text_repel(data = data_plot[!is.na(cross_hotspot_source) & 
        cross_hotspot_source == "K19", ], aes(x = distance_bp, y = mean, label = paste0(Pos_real, " (K19)")), nudge_y = 0.08, 
        size = text_repel_size * 0.8, fontface = "bold", color = "#C68EFD") + ggrepel::geom_text_repel(data = data_plot[!is.na(cross_hotspot_source) & 
        cross_hotspot_source == "RAF1", ], aes(x = distance_bp, y = mean, label = paste0(Pos_real, " (RAF1)")), nudge_y = 0.08, 
        size = text_repel_size * 0.8, fontface = "bold", color = "#09B636") + ggplot2::xlab(expression(paste("Distance to binding partner (" * 
        ring(A) * ")"))) + ggplot2::ylab("Weighted mean |ddG| (kcal/mol)") + ggplot2::labs(color = "Site Type", shape = "Secondary Structure") + 
        ggplot2::facet_wrap(~assay, ncol = 3) + ggplot2::scale_color_manual(values = color_values, labels = color_labels, 
        breaks = names(color_values), drop = FALSE) + ggplot2::scale_shape_manual(values = c(`beta strand` = 15, `alpha helix` = 16, 
        others = 17), drop = FALSE) + ggplot2::theme_classic(base_size = base_font_size) + ggplot2::theme(axis.text.x = element_text(size = base_font_size * 
        0.8, vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size = base_font_size * 0.8, vjust = 0.5, hjust = 0.5), 
        text = element_text(size = base_font_size), legend.position = "right", strip.text.x = element_text(size = base_font_size), 
        strip.background = element_rect(colour = "white", fill = "white"), panel.spacing = unit(0.2, "mm"), legend.text = element_text(size = base_font_size * 
            0.6), plot.margin = margin(0, 1, 0, 1, "mm"), legend.margin = margin(0, 0, 0, -2, "mm"), legend.spacing.y = unit(0, 
            "mm"), legend.key.height = unit(4, "mm"))
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, plot = p, device = cairo_pdf, height = 6, width = 16, dpi = 300)
        message(paste("Plot saved to:", output_file))
    }
    return(list(plot = p, data = data_plot, allosteric_list = allosteric_list, threshold = reg_threshold, cross_hotspots = data_plot[!is.na(cross_hotspot_source), 
        ]))
}

