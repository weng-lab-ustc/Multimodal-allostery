#' Plot Mapped Anticorrelated.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param k13_k19_data Value supplied for `k13_k19_data`.
#' @param raf1_k27_anticorrelated Value supplied for `raf1_k27_anticorrelated`.
#' @param nbp_residues Residue positions assigned to the nucleotide-binding pocket.
#' @param xlim Value supplied for `xlim`.
#' @param ylim Value supplied for `ylim`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_mapped_anticorrelated <- function(k13_k19_data, raf1_k27_anticorrelated, nbp_residues,
                                       xlim = c(-1.5, 3), ylim = c(-1.5, 3)) {
    merged_data <- copy(k13_k19_data$data)
    threshold_K13 <- k13_k19_data$threshold_K13
    threshold_K19 <- k13_k19_data$threshold_K19
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
    cor_test <- cor.test(merged_data$ddG_K13, merged_data$ddG_K19)
    r_value <- round(cor_test$estimate, 3)
    p_value <- cor_test$p.value
    sig_stars <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "")))
    r_label <- paste0("R = ", r_value, sig_stars)
    p <- ggplot() + theme_classic(base_size = 20) + geom_vline(xintercept = c(-threshold_K13, threshold_K13), linetype = "dashed", 
        color = "grey50") + geom_hline(yintercept = c(-threshold_K19, threshold_K19), linetype = "dashed", color = "grey50") + 
        geom_point(data = merged_data[is_target == FALSE], aes(ddG_K13, ddG_K19), color = "grey80", size = 2, alpha = 0.4, 
            shape = 16) + geom_point(data = target_non_nbp, aes(ddG_K13, ddG_K19), color = "#F1DD10", shape = 16, size = 3) + 
        geom_point(data = target_nbp, aes(ddG_K13, ddG_K19), color = "#F1DD10", shape = 17, size = 3.5) + annotate("text", 
        x = xlim[1], y = ylim[2], hjust = 0, vjust = 1, label = r_label, size = 6) + labs(x = expression("Binding" ~ Delta * 
        Delta * G ~ "(K13) (kcal/mol)"), y = expression("Binding" ~ Delta * Delta * G ~ "(K19) (kcal/mol)"), title = "RAF1 vs K27 anticorrelated mutations mapped onto K13 vs K19") + 
        coord_cartesian(xlim = xlim, ylim = ylim) + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(hjust = 0.5))
    return(p)
}

