#' plot mapped anticorrelated
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_mapped_anticorrelated()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param k13_k19_data Input or option used by this workflow; see Usage for the exact interface.
#' @param raf1_k27_anticorrelated Input or option used by this workflow; see Usage for the exact interface.
#' @param xlim Input or option used by this workflow; see Usage for the exact interface.
#' @param ylim Input or option used by this workflow; see Usage for the exact interface.
#' @param fixed_threshold Input or option used by this workflow; see Usage for the exact interface.
#' @param nbp_residues Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_mapped_anticorrelated <- function (k13_k19_data, raf1_k27_anticorrelated, xlim = c(-1.5, 3), ylim = c(-1.5, 3), fixed_threshold = 0.40000000000000002, 
    nbp_residues = c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)) 
{
    merged_data <- data.table::copy(k13_k19_data$data)
    threshold_K13 <- k13_k19_data$threshold_K13
    threshold_K19 <- k13_k19_data$threshold_K19
    if (is.na(threshold_K13) || is.null(threshold_K13)) {
        threshold_K13 <- fixed_threshold
        warning("threshold_K13 is NA/Null, using FIXED_THRESHOLD")
    }
    if (is.na(threshold_K19) || is.null(threshold_K19)) {
        threshold_K19 <- fixed_threshold
        warning("threshold_K19 is NA/Null, using FIXED_THRESHOLD")
    }
    cat("\nPlotting thresholds:\n")
    cat("  K13 threshold:", threshold_K13, "\n")
    cat("  K19 threshold:", threshold_K19, "\n")
    merged_data[, `:=`(mt, toupper(trimws(mt)))]
    raf1_k27_anticorrelated[, `:=`(mt, toupper(trimws(mt)))]
    merged_data[, `:=`(is_NBP, Pos_real %in% nbp_residues)]
    merged_data[, `:=`(is_target, mt %in% raf1_k27_anticorrelated$mt)]
    merged_data[, `:=`(target_NBP_status, FALSE)]
    for (i in 1:nrow(raf1_k27_anticorrelated)) {
        merged_data[mt == raf1_k27_anticorrelated$mt[i], `:=`(target_NBP_status, raf1_k27_anticorrelated$is_NBP[i])]
    }
    cat("\nMatched mutations:", merged_data[is_target == TRUE, .N], "\n")
    print(merged_data[is_target == TRUE, .(mt, Pos_real, is_NBP, target_NBP_status, ddG_K13, ddG_K19)])
    target_nbp <- merged_data[is_target == TRUE & target_NBP_status == TRUE]
    target_non_nbp <- merged_data[is_target == TRUE & target_NBP_status == FALSE]
    cor_test <- stats::cor.test(merged_data$ddG_K13, merged_data$ddG_K19)
    r_value <- round(cor_test$estimate, 3)
    p_value <- cor_test$p.value
    sig_stars <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.050000000000000003, "*", 
        "")))
    r_label <- paste0("R = ", r_value, sig_stars)
    p <- ggplot2::ggplot() + ggplot2::theme_classic(base_size = 20) + ggplot2::geom_vline(xintercept = c(-threshold_K13, 
        threshold_K13), linetype = "dashed", color = "grey50", size = 0.80000000000000004) + ggplot2::geom_hline(yintercept = c(-threshold_K19, 
        threshold_K19), linetype = "dashed", color = "grey50", size = 0.80000000000000004) + ggplot2::geom_point(data = merged_data[is_target == 
        FALSE], ggplot2::aes(ddG_K13, ddG_K19), color = "grey80", size = 2, alpha = 0.40000000000000002, shape = 16) + ggplot2::geom_point(data = target_non_nbp, 
        ggplot2::aes(ddG_K13, ddG_K19), color = "#F1DD10", shape = 16, size = 3) + ggplot2::geom_point(data = target_nbp, 
        ggplot2::aes(ddG_K13, ddG_K19), color = "#F1DD10", shape = 17, size = 3.5) + ggplot2::annotate("text", x = xlim[1] + 
        0.10000000000000001, y = ylim[2] - 0.10000000000000001, hjust = 0, vjust = 1, label = r_label, size = 6) + ggplot2::labs(x = expression("Binding" ~ 
        Delta * Delta * G ~ "(K13) (kcal/mol)"), y = expression("Binding" ~ Delta * Delta * G ~ "(K19) (kcal/mol)"), title = "RAF1 vs K27 anticorrelated mutations mapped onto K13 vs K19") + 
        ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", 
        fill = NA, linewidth = 0.80000000000000004), axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), 
        plot.title = ggplot2::element_text(hjust = 0.5))
    return(p)
}

