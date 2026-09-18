#' analyze correlation outliers
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_analyze_correlation_outliers()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param data Input or option used by this workflow; see Usage for the exact interface.
#' @param x_var Input or option used by this workflow; see Usage for the exact interface.
#' @param y_var Input or option used by this workflow; see Usage for the exact interface.
#' @param num_outliers Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_analyze_correlation_outliers <- function (data, x_var, y_var, num_outliers = 10) 
{
    complete_data <- data[stats::complete.cases(data[[x_var]], data[[y_var]]),]
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

