#' Analyze Correlation Outliers.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data Input data frame or data table.
#' @param x_var Value supplied for `x_var`.
#' @param y_var Value supplied for `y_var`.
#' @param num_outliers Value supplied for `num_outliers`.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_analyze_correlation_outliers <- function(data, x_var, y_var, num_outliers = 10) {
    complete_data <- data[stats::complete.cases(data[, c(x_var, y_var)]), ]
    formula <- stats::as.formula(paste(y_var, "~", x_var))
    fit <- stats::lm(formula, data = complete_data)
    complete_data$abs_residual <- abs(stats::resid(fit))
    complete_data <- complete_data %>% dplyr::mutate(outlier_rank = rank(-abs_residual, ties.method = "first"), is_outlier = outlier_rank <= 
        num_outliers)
    outliers <- complete_data %>% dplyr::filter(is_outlier == TRUE) %>% dplyr::arrange(dplyr::desc(abs_residual)) %>% dplyr::mutate(K13_error_range = 2 * 
        `std_kcal/mol.x`, K19_error_range = 2 * `std_kcal/mol.y`, likely_noise = abs_residual < (K13_error_range + K19_error_range)/2)
    return(list(data = complete_data, outliers = outliers, correlation = stats::cor.test(complete_data[[x_var]], complete_data[[y_var]], 
        use = "complete.obs"), fit = fit, num_outliers = num_outliers, noise_count = sum(outliers$likely_noise, na.rm = TRUE), 
        x_var = x_var, y_var = y_var))
}

