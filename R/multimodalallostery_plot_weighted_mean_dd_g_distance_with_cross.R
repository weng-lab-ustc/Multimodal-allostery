#' plot weighted mean dd g distance with cross
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_weighted_mean_dd_g_distance_with_cross()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param ddG_files Input or option used by this workflow; see Usage for the exact interface.
#' @param all_assays Input or option used by this workflow; see Usage for the exact interface.
#' @param plot_assays Input or option used by this workflow; see Usage for the exact interface.
#' @param anno_file Input or option used by this workflow; see Usage for the exact interface.
#' @param x_intercept Input or option used by this workflow; see Usage for the exact interface.
#' @param output_file Input or option used by this workflow; see Usage for the exact interface.
#' @param base_font_size Input or option used by this workflow; see Usage for the exact interface.
#' @param point_size Input or option used by this workflow; see Usage for the exact interface.
#' @param text_repel_size Input or option used by this workflow; see Usage for the exact interface.
#' @param gtp_assay_for_all Input or option used by this workflow; see Usage for the exact interface.
#' @param gtp_assay_for_plot Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_weighted_mean_dd_g_distance_with_cross <- function (ddG_files, all_assays, plot_assays, anno_file, x_intercept = 5, output_file = NULL, base_font_size = 12, point_size = 0.5, 
    text_repel_size = 3, gtp_assay_for_all = NULL, gtp_assay_for_plot = NULL) 
{
    assay_ligand <- c(RAF1 = "RAF1", K13 = "K13", K19 = "K19", K27 = "K27", K55 = "K55", RAL = "RALGDS", PI3 = "PI3KCG", SOS = "SOS1")
    rect_input <- data.frame(xstart = c(3, 15, 38, 51, 67, 77, 87, 109, 127, 139, 148), xend = c(9, 24, 44, 57, 73, 84, 104, 
        115, 136, 143, 166), col = c("b1", "a1", "b2", "b3", "a2", "b4", "a3", "b5", "a4", "b6", "a5"))
    anno <- .ma_read_table(anno_file)
    anno[, `:=`(Pos_real, Pos)]
    weighted_mean_ddG_all <- list()
    for (i in seq_along(all_assays)) {
        weighted_mean_ddG_all[[all_assays[i]]] <- krasddpcams::krasddpcams__get_weighted_mean_abs_ddG_mutcount(ddG = ddG_files[i], 
            assay_sele = all_assays[i])
    }
    weighted_mean_ddG_plot <- list()
    for (assayi in plot_assays) {
        idx <- which(all_assays == assayi)
        weighted_mean_ddG_plot[[assayi]] <- weighted_mean_ddG_all[[assayi]]
    }
    if (is.null(gtp_assay_for_all)) {
        gtp_assay_for_all <- all_assays[1]
    }
    if (is.null(gtp_assay_for_plot)) {
        gtp_assay_for_plot <- plot_assays[1]
    }
    cat(sprintf("Using '%s' for GTP binding site definition in ALL assays\n", gtp_assay_for_all))
    cat(sprintf("Using '%s' for GTP binding site definition in PLOT assays\n", gtp_assay_for_plot))
    data_all <- data.table::data.table()
    for (assayi in all_assays) {
        data_assayi <- merge(weighted_mean_ddG_all[[assayi]], anno, by = "Pos_real", all = TRUE)
        data_assayi[, `:=`(binding_type, "allosteric site")]
        data_assayi[get(paste0("scHAmin_ligand_", assay_ligand[[assayi]])) < x_intercept, `:=`(binding_type, "binding site")]
        data_assayi[, `:=`(binding_type_gtp_included, binding_type)]
        data_assayi[get(paste0("GXPMG_scHAmin_ligand_", assay_ligand[[gtp_assay_for_all]])) < x_intercept, `:=`(binding_type_gtp_included, 
            "GTP binding site")]
        data_all <- rbind(data_all, data_assayi)
    }
    reg_threshold <- data_all[binding_type == "binding site", sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = TRUE)/sum(1/.SD[[2]]^2, 
        na.rm = TRUE), .SDcols = c("mean", "sigma")]
    print(paste("Regression threshold (based on all", length(all_assays), "binders):", reg_threshold))
    data_plot <- data.table::data.table()
    for (assayi in plot_assays) {
        data_plot_assayi <- merge(weighted_mean_ddG_plot[[assayi]], anno, by = "Pos_real", all = TRUE)
        data_plot_assayi[, `:=`(binding_type, "allosteric site")]
        data_plot_assayi[get(paste0("scHAmin_ligand_", assay_ligand[[assayi]])) < x_intercept, `:=`(binding_type, "binding site")]
        data_plot_assayi[, `:=`(binding_type_gtp_included, binding_type)]
        data_plot_assayi[get(paste0("GXPMG_scHAmin_ligand_", gtp_assay_for_plot)) < x_intercept, `:=`(binding_type_gtp_included, 
            "GTP binding site")]
        data_plot <- rbind(data_plot, data_plot_assayi)
    }
    data_plot[, `:=`(site_type, "Reminder")]
    data_plot[binding_type_gtp_included == "binding site", `:=`(site_type, "Binding interface site")]
    data_plot[binding_type_gtp_included == "GTP binding site", `:=`(site_type, "Other GTP pocket site")]
    data_plot[binding_type_gtp_included == "GTP binding site" & mean > reg_threshold & binding_type != "binding site" & count > 
        9.5, `:=`(site_type, "Allosteric GTP pocket site")]
    data_plot[binding_type_gtp_included == "allosteric site" & mean > reg_threshold & count > 9.5, `:=`(site_type, "Major allosteric site")]
    data_plot[, `:=`(colors_type, "others")]
    rects_dt <- data.table::as.data.table(rect_input)
    for (b in c("b1", "b2", "b3", "b4", "b5", "b6")) {
        data_plot[Pos_real >= rects_dt[col == b, xstart] & Pos_real <= rects_dt[col == b, xend], `:=`(colors_type, b)]
    }
    for (a in c("a1", "a2", "a3", "a4", "a5")) {
        data_plot[Pos_real >= rects_dt[col == a, xstart] & Pos_real <= rects_dt[col == a, xend], `:=`(colors_type, a)]
    }
    data_plot[, `:=`(shape, "others")]
    data_plot[data.table::`%chin%`(colors_type, c("b1", "b2", "b3", "b4", "b5", "b6")), `:=`(shape, "beta strand")]
    data_plot[data.table::`%chin%`(colors_type, c("a1", "a2", "a3", "a4", "a5")), `:=`(shape, "alpha helix")]
    data_plot <- data_plot[Pos_real > 1 & count > 9.5, ]
    data_plot <- within(data_plot, site_type <- factor(site_type, levels = c("Binding interface site", "Allosteric GTP pocket site", 
        "Other GTP pocket site", "Major allosteric site", "Reminder")))
    data_plot <- within(data_plot, assay <- factor(assay, levels = plot_assays))
    allosteric_list <- list()
    for (assayi in plot_assays) {
        data_plot[assay == assayi, `:=`(distance_bp, get(paste0("scHAmin_ligand_", assay_ligand[[assayi]])))]
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
    color_values <- c(`Binding interface site` = "#F4270C", `Allosteric GTP pocket site` = "#1B38A6", `Major allosteric site` = "#F4AD0C", 
        `Other GTP pocket site` = "#75C2F6", Reminder = "gray")
    color_labels <- c(`Binding interface site` = "Binding interface", `Allosteric GTP pocket site` = "Allosteric GTP pocket", 
        `Major allosteric site` = "Major allosteric", `Other GTP pocket site` = "Other GTP pocket", Reminder = "Others")
    multimodalallostery_print_allosteric_statistics(data_plot = data_plot, allosteric_list = allosteric_list, assays = plot_assays, 
        reg_threshold = reg_threshold)
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = data_plot, mapping = ggplot2::aes(x = distance_bp, y = mean, color = site_type, 
        shape = as.factor(shape)), size = point_size) + ggplot2::geom_pointrange(data = data_plot, ggplot2::aes(x = distance_bp, 
        y = mean, color = site_type, ymin = mean - sigma, ymax = mean + sigma, shape = as.factor(shape)), size = point_size) + 
        ggplot2::geom_hline(yintercept = reg_threshold, linetype = 2, linewidth = 0.29999999999999999) + ggplot2::geom_vline(xintercept = x_intercept, 
        linetype = 2, linewidth = 0.29999999999999999) + ggplot2::geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.29999999999999999) + 
        ggplot2::geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.29999999999999999) + ggplot2::xlab(expression(paste("Distance to binding partner (" * 
        ring(A) * ")"))) + ggplot2::ylab("Weighted mean |ddG| (kcal/mol)") + ggplot2::labs(color = "Site Type", shape = "Secondary Structure") + 
        ggplot2::facet_wrap(~assay, ncol = 3) + ggplot2::scale_color_manual(values = color_values, labels = color_labels, 
        breaks = names(color_values), drop = FALSE) + ggplot2::scale_shape_manual(values = c(`beta strand` = 15, `alpha helix` = 16, 
        others = 17), drop = FALSE) + ggplot2::theme_classic(base_size = base_font_size) + ggplot2::theme(axis.text.x = ggplot2::element_text(size = base_font_size * 
        0.80000000000000004, vjust = 0.5, hjust = 0.5), axis.text.y = ggplot2::element_text(size = base_font_size * 0.80000000000000004, 
        vjust = 0.5, hjust = 0.5), text = ggplot2::element_text(size = base_font_size), legend.position = "right", strip.text.x = ggplot2::element_text(size = base_font_size), 
        strip.background = ggplot2::element_rect(colour = "white", fill = "white"), panel.spacing = grid::unit(0.20000000000000001, 
            "mm"), legend.text = ggplot2::element_text(size = base_font_size * 0.59999999999999998), plot.margin = ggplot2::margin(0, 
            1, 0, 1, "mm"), legend.margin = ggplot2::margin(0, 0, 0, -2, "mm"), legend.spacing.y = grid::unit(0, "mm"), legend.key.height = grid::unit(4, 
            "mm"))
    p <- p + ggrepel::geom_text_repel(data = data_plot[site_type == "Major allosteric site", ], ggplot2::aes(x = distance_bp, 
        y = mean, label = Pos_real), nudge_y = 0.050000000000000003, color = "#F4AD0C", size = text_repel_size, fontface = "bold")
    p <- p + ggrepel::geom_text_repel(data = data_plot[site_type == "Allosteric GTP pocket site", ], ggplot2::aes(x = distance_bp, 
        y = mean, label = Pos_real), nudge_y = 0.050000000000000003, color = "#1B38A6", size = text_repel_size, fontface = "bold")
    if (!is.null(output_file)) {
        output_dir <- dirname(output_file)
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        message("Saving plot using ggsave...")
        .ma_save_plot(output_file, plot = p, device = "pdf", height = 6, width = 16, dpi = 300)
        if (file.exists(output_file)) {
            file_size <- file.info(output_file)$size
            message(paste("Plot saved to:", output_file))
            message(paste("File size:", file_size, "bytes"))
        }
    }
    return(list(plot = p, data = data_plot, allosteric_list = allosteric_list, threshold = reg_threshold))
}

