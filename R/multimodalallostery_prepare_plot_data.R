#' Prepare Plot Data.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param NBP_results Value supplied for `NBP_results`.
#' @param second_shell_results Value supplied for `second_shell_results`.
#' @param beta_sheet_results Value supplied for `beta_sheet_results`.
#' @param assay_name Value supplied for `assay_name`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_plot_data <- function(NBP_results, second_shell_results, beta_sheet_results, assay_name) {
    total_allosteric_all_regions <- NBP_results$Case_in_NBP[1] + NBP_results$Case_out_NBP[1]
    nbp_data <- data.table(Region = "NBP", Category = c("Allosteric", "Inhibit", "Stabilize"), Cases_in_region = c(NBP_results$Case_in_NBP[1], 
        NBP_results$Case_in_NBP[2], NBP_results$Case_in_NBP[3]), OR = NBP_results$OR, P_value = NBP_results$P_value)
    nbp_data[, `:=`(Percentage, (Cases_in_region/total_allosteric_all_regions) * 100)]
    nbp_data[, `:=`(Error, {
        n <- total_allosteric_all_regions
        p <- Cases_in_region/n
        z <- 1.96
        lower <- (p + z^2/(2 * n) - z * sqrt((p * (1 - p)/n) + z^2/(4 * n^2)))/(1 + z^2/n)
        upper <- (p + z^2/(2 * n) + z * sqrt((p * (1 - p)/n) + z^2/(4 * n^2)))/(1 + z^2/n)
        (upper - lower)/2 * 100
    })]
    if (!is.null(second_shell_results) && nrow(second_shell_results) > 0) {
        shell_data <- data.table(Region = "Second shell", Category = c("Allosteric", "Inhibit", "Stabilize"), Cases_in_region = c(second_shell_results$Case_in_region[1], 
            second_shell_results$Case_in_region[2], second_shell_results$Case_in_region[3]), OR = second_shell_results$OR, 
            P_value = second_shell_results$P_value)
        shell_data[, `:=`(Percentage, (Cases_in_region/total_allosteric_all_regions) * 100)]
        shell_data[, `:=`(Error, {
            n <- total_allosteric_all_regions
            p <- Cases_in_region/n
            z <- 1.96
            lower <- (p + z^2/(2 * n) - z * sqrt((p * (1 - p)/n) + z^2/(4 * n^2)))/(1 + z^2/n)
            upper <- (p + z^2/(2 * n) + z * sqrt((p * (1 - p)/n) + z^2/(4 * n^2)))/(1 + z^2/n)
            (upper - lower)/2 * 100
        })]
    }
    else {
        shell_data <- NULL
    }
    if (!is.null(beta_sheet_results) && nrow(beta_sheet_results) > 0) {
        beta_data <- data.table(Region = "Beta sheet", Category = c("Allosteric", "Inhibit", "Stabilize"), Cases_in_region = c(beta_sheet_results$Case_in_region[1], 
            beta_sheet_results$Case_in_region[2], beta_sheet_results$Case_in_region[3]), OR = beta_sheet_results$OR, P_value = beta_sheet_results$P_value)
        beta_data[, `:=`(Percentage, (Cases_in_region/total_allosteric_all_regions) * 100)]
        beta_data[, `:=`(Error, {
            n <- total_allosteric_all_regions
            p <- Cases_in_region/n
            z <- 1.96
            lower <- (p + z^2/(2 * n) - z * sqrt((p * (1 - p)/n) + z^2/(4 * n^2)))/(1 + z^2/n)
            upper <- (p + z^2/(2 * n) + z * sqrt((p * (1 - p)/n) + z^2/(4 * n^2)))/(1 + z^2/n)
            (upper - lower)/2 * 100
        })]
    }
    else {
        beta_data <- NULL
    }
    plot_data <- rbind(nbp_data, shell_data, beta_data, fill = TRUE)
    plot_data[, `:=`(Region, factor(Region, levels = c("NBP", "Second shell", "Beta sheet")))]
    plot_data[, `:=`(Significance, ifelse(P_value < 0.05, ifelse(P_value < 0.01, ifelse(P_value < 0.001, "***", "**"), "*"), 
        " ns"))]
    plot_data[, `:=`(Sig_label, ifelse(OR > 1, paste0("OR = ", round(OR, 2), Significance), paste0("OR = ", round(OR, 2))))]
    return(plot_data)
}

