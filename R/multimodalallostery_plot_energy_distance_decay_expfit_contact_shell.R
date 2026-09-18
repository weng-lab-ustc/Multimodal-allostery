#' plot energy distance decay expfit contact shell
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_energy_distance_decay_expfit_contact_shell()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input Input or option used by this workflow; see Usage for the exact interface.
#' @param assay_sele Input or option used by this workflow; see Usage for the exact interface.
#' @param contact_shell Input or option used by this workflow; see Usage for the exact interface.
#' @param x_range Input or option used by this workflow; see Usage for the exact interface.
#' @param y_range Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_energy_distance_decay_expfit_contact_shell <- function (input, assay_sele, contact_shell, x_range = c(0, 10), y_range = c(0, 3)) 
{
    data <- .ma_read_table(input)
    data <- data[, `:=`(Pos_real, Pos + 1)]
    data <- data[, c(20:23)]
    colnames(data)[1:3] <- paste0(colnames(data)[1:3], "_", assay_sele)
    anno_final <- merge(contact_shell, data, by = "Pos_real", all = FALSE)
    x_col <- paste0(assay_sele, "_contact_shell")
    y_col <- paste0("mean_kcal/mol_", assay_sele)
    title <- paste0("Contact shell to ", assay_sele)
    xvector <- anno_final[[x_col]]
    yvector <- abs(anno_final[[y_col]])
    df <- data.frame(Pos_real = anno_final$Pos_real, x = xvector, y = yvector)
    df <- df[stats::complete.cases(df), ]
    df <- df[df$x > 0, ]
    df_per_residue <- df %>% dplyr::group_by(Pos_real, x) %>% dplyr::summarise(y_median = stats::median(y, na.rm = TRUE), 
        .groups = "drop")
    df_interface <- df_per_residue[df_per_residue$x == 1, ]
    df_shell <- df_per_residue[df_per_residue$x > 1, ]
    df_mut_fit <- df[df$x > 1, ]
    fit_mut <- tryCatch(stats::nls(y ~ a * exp(b * x), data = df_mut_fit, start = list(a = 1, b = -0.10000000000000001)), 
        error = function(e) NULL)
    fit_mut_df <- data.frame()
    annotation_mut <- NULL
    if (!is.null(fit_mut)) {
        x_seq <- seq(min(df_mut_fit$x), max(df_mut_fit$x), length.out = 200)
        fit_mut_df <- data.frame(x = x_seq, y = stats::predict(fit_mut, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_mut)$coefficients
        p_val_b <- coefs["b", "Pr(>|t|)"]
        annotation_mut <- paste0("Mutation-level:\n", "a = ", round(coefs["a", "Estimate"], 3), "\nb = ", round(coefs["b", 
            "Estimate"], 3), "\np = ", signif(p_val_b, 3))
    }
    df_res_fit <- df_shell[, c("x", "y_median")]
    colnames(df_res_fit)[2] <- "y"
    fit_res <- tryCatch(stats::nls(y ~ a * exp(b * x), data = df_res_fit, start = list(a = 1, b = -0.10000000000000001)), 
        error = function(e) NULL)
    fit_res_df <- data.frame()
    annotation_res <- NULL
    if (!is.null(fit_res)) {
        x_seq <- seq(min(df_res_fit$x), max(df_res_fit$x), length.out = 200)
        fit_res_df <- data.frame(x = x_seq, y = stats::predict(fit_res, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_res)$coefficients
        p_val_b <- coefs["b", "Pr(>|t|)"]
        annotation_res <- paste0("Residue-level:\n", "a = ", round(coefs["a", "Estimate"], 3), "\nb = ", round(coefs["b", 
            "Estimate"], 3), "\np = ", signif(p_val_b, 3))
    }
    df_median_shell <- df_shell
    colnames(df_median_shell)[3] <- "y"
    df_median_interface <- df_interface
    colnames(df_median_interface)[3] <- "y"
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = df[df$x > 1, ], ggplot2::aes(x = x, y = y), alpha = 0.10000000000000001, 
        size = 2, color = "#75C2F6") + ggplot2::geom_point(data = df[df$x == 1, ], ggplot2::aes(x = x, y = y), alpha = 0.10000000000000001, 
        size = 2, color = "#FFB6C1") + ggplot2::geom_point(data = df_median_shell, ggplot2::aes(x = x, y = y), color = "#1B38A6", 
        size = 2) + ggplot2::geom_point(data = df_median_interface, ggplot2::aes(x = x, y = y), color = "#8B0000", size = 2) + 
        ggplot2::geom_line(data = fit_mut_df, ggplot2::aes(x = x, y = y), color = "#75C2F6", linewidth = 1) + ggplot2::geom_line(data = fit_res_df, 
        ggplot2::aes(x = x, y = y), color = "#1B38A6", linewidth = 1, linetype = "dashed") + ggplot2::scale_x_continuous(limits = x_range, 
        expand = c(0, 0), breaks = seq(0, 10, by = 1)) + ggplot2::scale_y_continuous(limits = y_range, expand = c(0, 0), 
        breaks = seq(0, 3, by = 1)) + ggplot2::theme_classic(base_size = 10) + ggplot2::labs(x = title, y = paste0("Binding |\u0394\u0394G| (", 
        assay_sele, ") (kcal/mol)"))
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

