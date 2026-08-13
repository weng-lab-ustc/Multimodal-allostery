#' Plot Energy Distance Decay Expfit Contact Shell.
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
multimodalallostery_plot_energy_distance_decay_expfit_contact_shell <- function(input, assay_sele, contact_shell, x_range = c(0, 10), y_range = c(0, 
    3)) {
    data <- fread(input)
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
    df <- df[complete.cases(df), ]
    df <- df[df$x > 0, ]
    df_per_residue <- df %>% group_by(Pos_real, x) %>% summarise(y_median = median(y, na.rm = TRUE), .groups = "drop")
    df_interface <- df_per_residue[df_per_residue$x == 1, ]
    df_shell <- df_per_residue[df_per_residue$x > 1, ]
    df_mut_fit <- df[df$x > 1, ]
    fit_mut <- tryCatch(nls(y ~ a * exp(b * x), data = df_mut_fit, start = list(a = 1, b = -0.1)), error = function(e) NULL)
    fit_mut_df <- data.frame()
    annotation_mut <- NULL
    if (!is.null(fit_mut)) {
        x_seq <- seq(min(df_mut_fit$x), max(df_mut_fit$x), length.out = 200)
        fit_mut_df <- data.frame(x = x_seq, y = predict(fit_mut, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_mut)$coefficients
        p_val_b <- coefs["b", "Pr(>|t|)"]
        annotation_mut <- paste0("Mutation-level:\n", "a = ", round(coefs["a", "Estimate"], 3), "\nb = ", round(coefs["b", 
            "Estimate"], 3), "\np = ", signif(p_val_b, 3))
    }
    df_res_fit <- df_shell[, c("x", "y_median")]
    colnames(df_res_fit)[2] <- "y"
    fit_res <- tryCatch(nls(y ~ a * exp(b * x), data = df_res_fit, start = list(a = 1, b = -0.1)), error = function(e) NULL)
    fit_res_df <- data.frame()
    annotation_res <- NULL
    if (!is.null(fit_res)) {
        x_seq <- seq(min(df_res_fit$x), max(df_res_fit$x), length.out = 200)
        fit_res_df <- data.frame(x = x_seq, y = predict(fit_res, newdata = data.frame(x = x_seq)))
        coefs <- summary(fit_res)$coefficients
        p_val_b <- coefs["b", "Pr(>|t|)"]
        annotation_res <- paste0("Residue-level:\n", "a = ", round(coefs["a", "Estimate"], 3), "\nb = ", round(coefs["b", 
            "Estimate"], 3), "\np = ", signif(p_val_b, 3))
    }
    df_median_shell <- df_shell
    colnames(df_median_shell)[3] <- "y"
    df_median_interface <- df_interface
    colnames(df_median_interface)[3] <- "y"
    p <- ggplot() + geom_point(data = df[df$x > 1, ], aes(x = x, y = y), alpha = 0.1, size = 2, color = "#75C2F6") + geom_point(data = df[df$x == 
        1, ], aes(x = x, y = y), alpha = 0.1, size = 2, color = "#FFB6C1") + geom_point(data = df_median_shell, aes(x = x, 
        y = y), color = "#1B38A6", size = 2) + geom_point(data = df_median_interface, aes(x = x, y = y), color = "#8B0000", 
        size = 2) + geom_line(data = fit_mut_df, aes(x = x, y = y), color = "#75C2F6", linewidth = 1) + geom_line(data = fit_res_df, 
        aes(x = x, y = y), color = "#1B38A6", linewidth = 1, linetype = "dashed") + scale_x_continuous(limits = x_range, 
        expand = c(0, 0), breaks = seq(0, 10, by = 1)) + scale_y_continuous(limits = y_range, expand = c(0, 0), breaks = seq(0, 
        3, by = 1)) + theme_classic(base_size = 10) + labs(x = title, y = paste0("Binding |\u0394\u0394G| (", assay_sele, ") (kcal/mol)"))
    if (!is.null(annotation_mut)) {
        p <- p + annotate("text", x = max(x_range) * 0.95, y = max(y_range) * 0.95, label = annotation_mut, hjust = 1, vjust = 1, 
            size = 2.8, color = "#75C2F6")
    }
    if (!is.null(annotation_res)) {
        p <- p + annotate("text", x = max(x_range) * 0.95, y = max(y_range) * 0.6, label = annotation_res, hjust = 1, vjust = 1, 
            size = 2.8, color = "#1B38A6")
    }
    return(p)
}

