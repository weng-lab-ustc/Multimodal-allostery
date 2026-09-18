#' plot energy distance decay expfit annotation
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_energy_distance_decay_expfit()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_sele Input or option used by this workflow; see Usage for the exact interface.
#' @param anno_file Input or option used by this workflow; see Usage for the exact interface.
#' @param x_range Input or option used by this workflow; see Usage for the exact interface.
#' @param y_range Input or option used by this workflow; see Usage for the exact interface.
#' @param plot_width Input or option used by this workflow; see Usage for the exact interface.
#' @param plot_height Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_energy_distance_decay_expfit <- function (input, assay_sele, anno_file, x_range = c(0, 35), y_range = c(0, 3), plot_width = 4, plot_height = 4) 
{
    data <- .ma_read_table(input)
    data <- data[, `:=`(Pos_real, Pos + 1)]
    data <- data[, c(20:23)]
    colnames(data)[1:3] <- paste0(colnames(data)[1:3], "_", assay_sele)
    anno <- .ma_read_table(anno_file)
    anno <- anno[, `:=`(Pos_real, Pos)]
    anno_final <- merge(anno, data, by = "Pos_real", all = FALSE)
    x_col <- paste0("scHAmin_ligand_", assay_sele)
    y_col <- paste0("mean_kcal/mol_", assay_sele)
    title <- paste0("Distance to ", assay_sele, " (\u00c5)")
    xvector <- anno_final[[x_col]]
    yvector <- abs(anno_final[[y_col]])
    df <- data.frame(x = xvector, y = yvector)
    df <- df[stats::complete.cases(df), ]
    df <- df[df$x > 0, ]
    df_bi <- df[df$x < 5, ]
    df_non_bi <- df[df$x >= 5, ]
    df_residue <- df %>% dplyr::group_by(x) %>% dplyr::summarise(y = stats::median(y, na.rm = TRUE), .groups = "drop")
    df_residue_non_bi <- df_residue[df_residue$x >= 5, ]
    fit_mut <- tryCatch(stats::nls(y ~ a * exp(b * x), data = df_non_bi, start = list(a = 1, b = -0.10000000000000001)), 
        error = function(e) NULL)
    fit_mut_df <- data.frame()
    annotation_mut <- NULL
    if (!is.null(fit_mut)) {
        x_seq <- seq(min(df_non_bi$x), max(df_non_bi$x), length.out = 200)
        fit_mut_df <- data.frame(x = x_seq, y = stats::predict(fit_mut, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_mut)$coefficients
        p_val_b <- coefs["b", "Pr(>|t|)"]
        annotation_mut <- paste0("Mutation-level:\n", "a = ", round(coefs["a", "Estimate"], 3), "\nb = ", round(coefs["b", 
            "Estimate"], 3), "\np = ", signif(p_val_b, 3))
    }
    fit_res <- tryCatch(stats::nls(y ~ a * exp(b * x), data = df_residue_non_bi, start = list(a = 1, b = -0.10000000000000001)), 
        error = function(e) NULL)
    fit_res_df <- data.frame()
    annotation_res <- NULL
    if (!is.null(fit_res)) {
        x_seq <- seq(min(df_residue_non_bi$x), max(df_residue_non_bi$x), length.out = 200)
        fit_res_df <- data.frame(x = x_seq, y = stats::predict(fit_res, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_res)$coefficients
        p_val_b <- coefs["b", "Pr(>|t|)"]
        annotation_res <- paste0("Residue-level:\n", "a = ", round(coefs["a", "Estimate"], 3), "\nb = ", round(coefs["b", 
            "Estimate"], 3), "\np = ", signif(p_val_b, 3))
    }
    df_median_non_bi <- df_non_bi %>% dplyr::group_by(x) %>% dplyr::summarise(y = stats::median(y, na.rm = TRUE), .groups = "drop")
    df_median_bi <- df_bi %>% dplyr::group_by(x) %>% dplyr::summarise(y = stats::median(y, na.rm = TRUE), .groups = "drop")
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = df_non_bi, ggplot2::aes(x = x, y = y), alpha = 0.10000000000000001, 
        size = 1.5, color = "#75C2F6") + ggplot2::geom_point(data = df_bi, ggplot2::aes(x = x, y = y), alpha = 0.10000000000000001, 
        size = 1.5, color = "#FFB6C1") + ggplot2::geom_point(data = df_median_non_bi, ggplot2::aes(x = x, y = y), color = "#1B38A6", 
        size = 2) + ggplot2::geom_point(data = df_median_bi, ggplot2::aes(x = x, y = y), color = "#8B0000", size = 2) + ggplot2::geom_line(data = fit_mut_df, 
        ggplot2::aes(x = x, y = y), color = "#75C2F6", linewidth = 0.80000000000000004) + ggplot2::geom_line(data = fit_res_df, 
        ggplot2::aes(x = x, y = y), color = "#1B38A6", linewidth = 0.80000000000000004, linetype = "dashed") + ggplot2::geom_vline(xintercept = 5, 
        linetype = "dashed", color = "gray50", linewidth = 0.5) + ggplot2::scale_x_continuous(limits = x_range, expand = c(0, 
        0)) + ggplot2::scale_y_continuous(limits = y_range, expand = c(0, 0)) + ggplot2::theme_classic(base_size = 10) + 
        ggplot2::labs(x = title, y = paste0("Binding |\u0394\u0394G| (", assay_sele, ") (kcal/mol)")) + ggplot2::theme(axis.title = ggplot2::element_text(size = 10), 
        axis.text = ggplot2::element_text(size = 10), plot.title = ggplot2::element_text(size = 10))
    if (!is.null(annotation_mut)) {
        p <- p + ggplot2::annotate("text", x = max(x_range) * 0.94999999999999996, y = max(y_range) * 0.94999999999999996, 
            label = annotation_mut, hjust = 1, vjust = 1, size = 2.7999999999999998, color = "#75C2F6")
    }
    if (!is.null(annotation_res)) {
        p <- p + ggplot2::annotate("text", x = max(x_range) * 0.94999999999999996, y = max(y_range) * 0.59999999999999998, 
            label = annotation_res, hjust = 1, vjust = 1, size = 2.7999999999999998, color = "#1B38A6")
    }
    return(p)
}

