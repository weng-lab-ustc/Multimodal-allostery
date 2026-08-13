#' Plot Energy Distance Decay Expfit Directional No Fdr.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input Input file path or input data object.
#' @param assay_sele Selected assay identifier.
#' @param contact_shell Value supplied for `contact_shell`.
#' @param x_range Value supplied for `x_range`.
#' @param y_range Value supplied for `y_range`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_energy_distance_decay_expfit_directional_no_fdr <- function(input, assay_sele, contact_shell, x_range = c(0, 10), y_range = c(-1.5, 
    3)) {
    data <- fread(input)
    data <- data[, `:=`(Pos_real, Pos + 1)]
    data <- data[, c(20:23)]
    colnames(data)[1:3] <- paste0(colnames(data)[1:3], "_", assay_sele)
    anno_final <- merge(contact_shell, data, by = "Pos_real", all = FALSE)
    x_col <- paste0(assay_sele, "_contact_shell")
    mean_col <- paste0("mean_kcal/mol_", assay_sele)
    df_inhibit <- anno_final[get(mean_col) > 0, .(x = get(x_col), y = abs(get(mean_col)))]
    df_stabilize <- anno_final[get(mean_col) < 0, .(x = get(x_col), y = get(mean_col))]
    df_inhibit <- df_inhibit[complete.cases(df_inhibit) & is.finite(x) & is.finite(y) & x > 0]
    df_stabilize <- df_stabilize[complete.cases(df_stabilize) & is.finite(x) & is.finite(y) & x > 0]
    df_inhibit_fit <- df_inhibit[x > 1]
    df_stabilize_fit <- df_stabilize[x > 1]
    fit_inhibit <- tryCatch(nls(y ~ a * exp(b * x), data = df_inhibit_fit, start = list(a = 1, b = -0.1)), error = function(e) NULL)
    fit_inhibit_df <- data.frame()
    annotation_inhibit <- NULL
    if (!is.null(fit_inhibit)) {
        x_seq <- seq(min(df_inhibit_fit$x), max(df_inhibit_fit$x), length.out = 200)
        fit_inhibit_df <- data.frame(x = x_seq, y = predict(fit_inhibit, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_inhibit)$coefficients
        annotation_inhibit <- paste0("Destabilize (\u0394\u0394G>0):\n", "a = ", round(coefs["a", "Estimate"], 3), "\nb = ", round(coefs["b", 
            "Estimate"], 3))
    }
    fit_stabilize <- tryCatch(nls(y ~ a * exp(b * x), data = df_stabilize_fit, start = list(a = -1, b = -0.1)), error = function(e) NULL)
    fit_stabilize_df <- data.frame()
    annotation_stabilize <- NULL
    if (!is.null(fit_stabilize)) {
        x_seq <- seq(min(df_stabilize_fit$x), max(df_stabilize_fit$x), length.out = 200)
        fit_stabilize_df <- data.frame(x = x_seq, y = predict(fit_stabilize, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_stabilize)$coefficients
        annotation_stabilize <- paste0("Stabilize (\u0394\u0394G<0):\n", "a = ", round(coefs["a", "Estimate"], 3), "\nb = ", round(coefs["b", 
            "Estimate"], 3))
    }
    df_inhibit_median <- df_inhibit %>% group_by(x) %>% summarise(y_median = median(y, na.rm = TRUE), .groups = "drop")
    df_stabilize_median <- df_stabilize %>% group_by(x) %>% summarise(y_median = median(y, na.rm = TRUE), .groups = "drop")
    p <- ggplot() + geom_point(data = df_inhibit, aes(x = x, y = y), alpha = 0.15, size = 1.5, color = "#F4270C") + geom_point(data = df_stabilize, 
        aes(x = x, y = y), alpha = 0.15, size = 1.5, color = "#1B38A6") + geom_point(data = df_inhibit_median, aes(x = x, 
        y = y_median), color = "#F4270C", size = 2, shape = 16) + geom_point(data = df_stabilize_median, aes(x = x, y = y_median), 
        color = "#1B38A6", size = 2, shape = 16) + geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) + 
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.3) + geom_line(data = fit_inhibit_df, 
        aes(x = x, y = y), color = "gray40", linewidth = 1) + geom_line(data = fit_stabilize_df, aes(x = x, y = y), color = "gray40", 
        linewidth = 1) + scale_x_continuous(limits = x_range, expand = c(0, 0), breaks = seq(0, 10, by = 2)) + scale_y_continuous(limits = y_range, 
        expand = c(0, 0), breaks = seq(-3, 3, by = 1)) + theme_classic(base_size = 10) + labs(x = paste0("Contact shell distance to ", 
        assay_sele, " (\u00C5)"), y = paste0("Binding \u0394\u0394G (", assay_sele, ") (kcal/mol)")) + theme(axis.title = element_text(size = 10), 
        axis.text = element_text(size = 10))
    if (!is.null(annotation_inhibit)) {
        p <- p + annotate("text", x = max(x_range) * 0.95, y = max(y_range) * 0.95, label = annotation_inhibit, hjust = 1, 
            vjust = 1, size = 2.5, color = "#F4270C")
    }
    if (!is.null(annotation_stabilize)) {
        p <- p + annotate("text", x = max(x_range) * 0.95, y = max(y_range) * 0.65, label = annotation_stabilize, hjust = 1, 
            vjust = 1, size = 2.5, color = "#1B38A6")
    }
    return(p)
}

