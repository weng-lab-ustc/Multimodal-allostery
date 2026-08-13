#' Plot Mapped both Promoting K27 K13.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param k27_k13_data Value supplied for `k27_k13_data`.
#' @param both_promoting_muts Value supplied for `both_promoting_muts`.
#' @param xlim Value supplied for `xlim`.
#' @param ylim Value supplied for `ylim`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_mapped_both_promoting_k27_k13 <- function(k27_k13_data, both_promoting_muts, xlim = c(-1.5, 3), ylim = c(-1.5, 3)) {
    merged_data <- copy(k27_k13_data$data)
    threshold_K27 <- k27_k13_data$threshold_K27
    threshold_K13 <- k27_k13_data$threshold_K13
    merged_data[, `:=`(mt, toupper(trimws(mt)))]
    both_promoting_muts[, `:=`(mt, toupper(trimws(mt)))]
    merged_data[, `:=`(is_target, mt %in% both_promoting_muts$mt)]
    merged_data[, `:=`(target_NBP_status, FALSE)]
    for (i in 1:nrow(both_promoting_muts)) {
        merged_data[mt == both_promoting_muts$mt[i], `:=`(target_NBP_status, both_promoting_muts$is_NBP[i])]
    }
    cat("\nMapped to K27 vs K13 - Matched mutations:", merged_data[is_target == TRUE, .N], "\n")
    if (merged_data[is_target == TRUE, .N] > 0) {
        print(merged_data[is_target == TRUE, .(mt, target_NBP_status, ddG_K27, ddG_K13)])
    }
    target_nbp <- merged_data[is_target == TRUE & target_NBP_status == TRUE]
    target_non_nbp <- merged_data[is_target == TRUE & target_NBP_status == FALSE]
    cor_test <- cor.test(merged_data$ddG_K27, merged_data$ddG_K13)
    r_value <- round(cor_test$estimate, 3)
    p_value <- cor_test$p.value
    sig_stars <- ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "")))
    r_label <- paste0("R = ", r_value, sig_stars)
    p <- ggplot() + theme_classic(base_size = 20) + geom_vline(xintercept = c(-threshold_K27, threshold_K27), linetype = "dashed", 
        color = "grey50") + geom_hline(yintercept = c(-threshold_K13, threshold_K13), linetype = "dashed", color = "grey50") + 
        geom_point(data = merged_data[is_target == FALSE], aes(ddG_K27, ddG_K13), color = "grey80", size = 2, alpha = 0.4, 
            shape = 16) + geom_point(data = target_non_nbp, aes(ddG_K27, ddG_K13), color = "#FFB0A5", shape = 16, size = 3) + 
        geom_point(data = target_nbp, aes(ddG_K27, ddG_K13), color = "#FFB0A5", shape = 17, size = 3.5) + annotate("text", 
        x = xlim[1], y = ylim[2], hjust = 0, vjust = 1, label = r_label, size = 6) + labs(x = expression("Binding" ~ Delta * 
        Delta * G ~ "(K27) (kcal/mol)"), y = expression("Binding" ~ Delta * Delta * G ~ "(K13) (kcal/mol)"), title = "K13 vs K19 both promoting mutations mapped onto K27 vs K13") + 
        coord_cartesian(xlim = xlim, ylim = ylim) + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(hjust = 0.5))
    return(p)
}

