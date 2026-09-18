#' plot energy distance decay expfit directional no fdr
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_energy_distance_decay_expfit_contact_shell_directional()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_sele Input or option used by this workflow; see Usage for the exact interface.
#' @param contact_shell Input or option used by this workflow; see Usage for the exact interface.
#' @param x_range Input or option used by this workflow; see Usage for the exact interface.
#' @param y_range Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_energy_distance_decay_expfit_contact_shell_directional <- function (input, assay_sele, contact_shell, x_range = c(0, 10), y_range = c(-1.5, 3)) 
{
    data <- .ma_read_table(input)
    data <- data[, `:=`(Pos_real, Pos + 1)]
    data <- data[, c(20:23)]
    colnames(data)[1:3] <- paste0(colnames(data)[1:3], "_", assay_sele)
    anno_final <- merge(contact_shell, data, by = "Pos_real", all = FALSE)
    x_col <- paste0(assay_sele, "_contact_shell")
    mean_col <- paste0("mean_kcal/mol_", assay_sele)
    df_inhibit <- anno_final[get(mean_col) > 0, .(x = get(x_col), y = abs(get(mean_col)))]
    df_stabilize <- anno_final[get(mean_col) < 0, .(x = get(x_col), y = get(mean_col))]
    df_inhibit <- df_inhibit[stats::complete.cases(df_inhibit) & is.finite(x) & is.finite(y) & x > 0]
    df_stabilize <- df_stabilize[stats::complete.cases(df_stabilize) & is.finite(x) & is.finite(y) & x > 0]
    df_inhibit_fit <- df_inhibit[x > 1]
    df_stabilize_fit <- df_stabilize[x > 1]
    fit_inhibit <- tryCatch(stats::nls(y ~ a * exp(b * x), data = df_inhibit_fit, start = list(a = 1, b = -0.10000000000000001)), 
        error = function(e) NULL)
    fit_inhibit_df <- data.frame()
    annotation_inhibit <- NULL
    if (!is.null(fit_inhibit)) {
        x_seq <- seq(min(df_inhibit_fit$x), max(df_inhibit_fit$x), length.out = 200)
        fit_inhibit_df <- data.frame(x = x_seq, y = stats::predict(fit_inhibit, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_inhibit)$coefficients
        annotation_inhibit <- paste0("Destabilize (\u0394\u0394G>0):\n", "a = ", round(coefs["a", "Estimate"], 3), "\nb = ", round(coefs["b", 
            "Estimate"], 3))
    }
    fit_stabilize <- tryCatch(stats::nls(y ~ a * exp(b * x), data = df_stabilize_fit, start = list(a = -1, b = -0.10000000000000001)), 
        error = function(e) NULL)
    fit_stabilize_df <- data.frame()
    annotation_stabilize <- NULL
    if (!is.null(fit_stabilize)) {
        x_seq <- seq(min(df_stabilize_fit$x), max(df_stabilize_fit$x), length.out = 200)
        fit_stabilize_df <- data.frame(x = x_seq, y = stats::predict(fit_stabilize, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_stabilize)$coefficients
        annotation_stabilize <- paste0("Stabilize (\u0394\u0394G<0):\n", "a = ", round(coefs["a", "Estimate"], 3), "\nb = ", round(coefs["b", 
            "Estimate"], 3))
    }
    df_inhibit_median <- df_inhibit %>% dplyr::group_by(x) %>% dplyr::summarise(y_median = stats::median(y, na.rm = TRUE), 
        .groups = "drop")
    df_stabilize_median <- df_stabilize %>% dplyr::group_by(x) %>% dplyr::summarise(y_median = stats::median(y, na.rm = TRUE), 
        .groups = "drop")
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = df_inhibit, ggplot2::aes(x = x, y = y), alpha = 0.14999999999999999, 
        size = 1.5, color = "#F4270C") + ggplot2::geom_point(data = df_stabilize, ggplot2::aes(x = x, y = y), alpha = 0.14999999999999999, 
        size = 1.5, color = "#1B38A6") + ggplot2::geom_point(data = df_inhibit_median, ggplot2::aes(x = x, y = y_median), 
        color = "#F4270C", size = 2, shape = 16) + ggplot2::geom_point(data = df_stabilize_median, ggplot2::aes(x = x, y = y_median), 
        color = "#1B38A6", size = 2, shape = 16) + ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", 
        linewidth = 0.5) + ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.29999999999999999) + 
        ggplot2::geom_line(data = fit_inhibit_df, ggplot2::aes(x = x, y = y), color = "gray40", linewidth = 1) + ggplot2::geom_line(data = fit_stabilize_df, 
        ggplot2::aes(x = x, y = y), color = "gray40", linewidth = 1) + ggplot2::scale_x_continuous(limits = x_range, expand = c(0, 
        0), breaks = seq(0, 10, by = 2)) + ggplot2::scale_y_continuous(limits = y_range, expand = c(0, 0), breaks = seq(-3, 
        3, by = 1)) + ggplot2::theme_classic(base_size = 10) + ggplot2::labs(x = paste0("Contact shell distance to ", assay_sele, 
        " (\u00c5)"), y = paste0("Binding \u0394\u0394G (", assay_sele, ") (kcal/mol)")) + ggplot2::theme(axis.title = ggplot2::element_text(size = 10), 
        axis.text = ggplot2::element_text(size = 10))
    if (!is.null(annotation_inhibit)) {
        p <- p + ggplot2::annotate("text", x = max(x_range) * 0.94999999999999996, y = max(y_range) * 0.94999999999999996, 
            label = annotation_inhibit, hjust = 1, vjust = 1, size = 2.5, color = "#F4270C")
    }
    if (!is.null(annotation_stabilize)) {
        p <- p + ggplot2::annotate("text", x = max(x_range) * 0.94999999999999996, y = max(y_range) * 0.65000000000000002, 
            label = annotation_stabilize, hjust = 1, vjust = 1, size = 2.5, color = "#1B38A6")
    }
    return(p)
}

