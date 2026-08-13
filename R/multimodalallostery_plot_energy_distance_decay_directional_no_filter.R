#' Plot Energy Distance Decay Directional No Filter.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input Input file path or input data object.
#' @param assay_sele Selected assay identifier.
#' @param anno_file Value supplied for `anno_file`.
#' @param x_range Value supplied for `x_range`.
#' @param y_range Value supplied for `y_range`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_energy_distance_decay_directional_no_filter <- function(input, assay_sele, anno_file, x_range = c(0, 35), y_range = c(-1.5, 
    3)) {
    data <- fread(input)
    data[, `:=`(Pos_real, Pos + 1)]
    anno <- fread(anno_file)
    anno[, `:=`(Pos_real, Pos)]
    anno_final <- merge(anno, data, by = "Pos_real", all = FALSE)
    x_col <- paste0("scHAmin_ligand_", assay_sele)
    df_inhibit <- anno_final[`mean_kcal/mol` > 0, .(x = get(x_col), y = abs(`mean_kcal/mol`))]
    df_activate <- anno_final[`mean_kcal/mol` < 0, .(x = get(x_col), y = `mean_kcal/mol`)]
    df_inhibit <- df_inhibit[complete.cases(df_inhibit)]
    df_activate <- df_activate[complete.cases(df_activate)]
    df_inhibit_fit <- df_inhibit[x >= 5]
    df_activate_fit <- df_activate[x >= 5]
    fit_inhibit <- tryCatch(nls(y ~ a * exp(b * x), data = df_inhibit_fit, start = list(a = 1, b = -0.1)), error = function(e) NULL)
    fit_inhibit_df <- data.frame()
    annotation_inhibit <- NULL
    if (!is.null(fit_inhibit)) {
        x_seq <- seq(min(df_inhibit_fit$x), max(df_inhibit_fit$x), length.out = 200)
        fit_inhibit_df <- data.frame(x = x_seq, y = predict(fit_inhibit, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_inhibit)$coefficients
        annotation_inhibit <- paste0("Inhibit binding\n", "a = ", round(coefs["a", "Estimate"], 3), "\n", "b = ", round(coefs["b", 
            "Estimate"], 3))
    }
    fit_activate <- tryCatch(nls(y ~ a * exp(b * x), data = df_activate_fit, start = list(a = -1, b = -0.1)), error = function(e) NULL)
    fit_activate_df <- data.frame()
    annotation_activate <- NULL
    if (!is.null(fit_activate)) {
        x_seq <- seq(min(df_activate_fit$x), max(df_activate_fit$x), length.out = 200)
        fit_activate_df <- data.frame(x = x_seq, y = predict(fit_activate, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_activate)$coefficients
        annotation_activate <- paste0("Stabilize binding\n", "a = ", round(coefs["a", "Estimate"], 3), "\n", "b = ", round(coefs["b", 
            "Estimate"], 3))
    }
    p <- ggplot() + geom_point(data = df_inhibit, aes(x = x, y = y), color = "#F4270C", alpha = 0.25, size = 1.5) + geom_point(data = df_activate, 
        aes(x = x, y = y), color = "#1B38A6", alpha = 0.25, size = 1.5) + geom_line(data = fit_inhibit_df, aes(x = x, y = y), 
        color = "gray40", linewidth = 1) + geom_line(data = fit_activate_df, aes(x = x, y = y), color = "gray40", linewidth = 1) + 
        geom_vline(xintercept = 5, linetype = "dashed", color = "gray50") + scale_x_continuous(limits = x_range, expand = c(0, 
        0)) + scale_y_continuous(limits = y_range, expand = c(0, 0)) + theme_classic(base_size = 10) + labs(x = paste0("Distance to ", 
        assay_sele, " (\u00C5)"), y = paste0("Binding \u0394\u0394G (", assay_sele, ") (kcal/mol)")) + theme(axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10))
    if (!is.null(annotation_inhibit)) {
        p <- p + annotate("text", x = max(x_range) * 0.95, y = max(y_range) * 0.95, label = annotation_inhibit, hjust = 1, 
            vjust = 1, size = 2.8, color = "#F4270C")
    }
    if (!is.null(annotation_activate)) {
        p <- p + annotate("text", x = max(x_range) * 0.95, y = max(y_range) * 0.65, label = annotation_activate, hjust = 1, 
            vjust = 1, size = 2.8, color = "#1B38A6")
    }
    return(p)
}

