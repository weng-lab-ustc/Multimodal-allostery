#' plot energy distance decay directional no filter
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_energy_distance_decay_expfit_directional()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_sele Input or option used by this workflow; see Usage for the exact interface.
#' @param anno_file Input or option used by this workflow; see Usage for the exact interface.
#' @param x_range Input or option used by this workflow; see Usage for the exact interface.
#' @param y_range Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_energy_distance_decay_expfit_directional <- function (input, assay_sele, anno_file, x_range = c(0, 35), y_range = c(-1.5, 3)) 
{
    data <- .ma_read_table(input)
    data[, `:=`(Pos_real, Pos + 1)]
    anno <- .ma_read_table(anno_file)
    anno[, `:=`(Pos_real, Pos)]
    anno_final <- merge(anno, data, by = "Pos_real", all = FALSE)
    x_col <- paste0("scHAmin_ligand_", assay_sele)
    df_inhibit <- anno_final[`mean_kcal/mol` > 0, .(x = get(x_col), y = abs(`mean_kcal/mol`))]
    df_activate <- anno_final[`mean_kcal/mol` < 0, .(x = get(x_col), y = `mean_kcal/mol`)]
    df_inhibit <- df_inhibit[stats::complete.cases(df_inhibit)]
    df_activate <- df_activate[stats::complete.cases(df_activate)]
    df_inhibit_fit <- df_inhibit[x >= 5]
    df_activate_fit <- df_activate[x >= 5]
    fit_inhibit <- tryCatch(stats::nls(y ~ a * exp(b * x), data = df_inhibit_fit, start = list(a = 1, b = -0.10000000000000001)), 
        error = function(e) NULL)
    fit_inhibit_df <- data.frame()
    annotation_inhibit <- NULL
    if (!is.null(fit_inhibit)) {
        x_seq <- seq(min(df_inhibit_fit$x), max(df_inhibit_fit$x), length.out = 200)
        fit_inhibit_df <- data.frame(x = x_seq, y = stats::predict(fit_inhibit, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_inhibit)$coefficients
        annotation_inhibit <- paste0("Inhibit binding\n", "a = ", round(coefs["a", "Estimate"], 3), "\n", "b = ", round(coefs["b", 
            "Estimate"], 3))
    }
    fit_activate <- tryCatch(stats::nls(y ~ a * exp(b * x), data = df_activate_fit, start = list(a = -1, b = -0.10000000000000001)), 
        error = function(e) NULL)
    fit_activate_df <- data.frame()
    annotation_activate <- NULL
    if (!is.null(fit_activate)) {
        x_seq <- seq(min(df_activate_fit$x), max(df_activate_fit$x), length.out = 200)
        fit_activate_df <- data.frame(x = x_seq, y = stats::predict(fit_activate, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_activate)$coefficients
        annotation_activate <- paste0("Stabilize binding\n", "a = ", round(coefs["a", "Estimate"], 3), "\n", "b = ", round(coefs["b", 
            "Estimate"], 3))
    }
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = df_inhibit, ggplot2::aes(x = x, y = y), color = "#F4270C", alpha = 0.25, 
        size = 1.5) + ggplot2::geom_point(data = df_activate, ggplot2::aes(x = x, y = y), color = "#1B38A6", alpha = 0.25, 
        size = 1.5) + ggplot2::geom_line(data = fit_inhibit_df, ggplot2::aes(x = x, y = y), color = "gray40", linewidth = 1) + 
        ggplot2::geom_line(data = fit_activate_df, ggplot2::aes(x = x, y = y), color = "gray40", linewidth = 1) + ggplot2::geom_vline(xintercept = 5, 
        linetype = "dashed", color = "gray50") + ggplot2::scale_x_continuous(limits = x_range, expand = c(0, 0)) + ggplot2::scale_y_continuous(limits = y_range, 
        expand = c(0, 0)) + ggplot2::theme_classic(base_size = 10) + ggplot2::labs(x = paste0("Distance to ", assay_sele, 
        " (\u00c5)"), y = paste0("Binding \u0394\u0394G (", assay_sele, ") (kcal/mol)")) + ggplot2::theme(axis.title = ggplot2::element_text(size = 10), 
        axis.text = ggplot2::element_text(size = 10))
    if (!is.null(annotation_inhibit)) {
        p <- p + ggplot2::annotate("text", x = max(x_range) * 0.94999999999999996, y = max(y_range) * 0.94999999999999996, 
            label = annotation_inhibit, hjust = 1, vjust = 1, size = 2.7999999999999998, color = "#F4270C")
    }
    if (!is.null(annotation_activate)) {
        p <- p + ggplot2::annotate("text", x = max(x_range) * 0.94999999999999996, y = max(y_range) * 0.65000000000000002, 
            label = annotation_activate, hjust = 1, vjust = 1, size = 2.7999999999999998, color = "#1B38A6")
    }
    return(p)
}

