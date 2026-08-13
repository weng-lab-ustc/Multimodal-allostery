#' Calculate Correlation Matrix.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param input_files Value supplied for `input_files`.
#' @param binder_order Value supplied for `binder_order`.
#' @param binding_sites_map Named assay-to-binding-site residue mapping.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_calculate_correlation_matrix <- function(input_files, binder_order, binding_sites_map) {
    assays <- names(input_files)
    n <- length(assays)
    assays <- binder_order[binder_order %in% assays]
    n <- length(assays)
    cor_matrix <- matrix(NA, nrow = n, ncol = n)
    rownames(cor_matrix) <- assays
    colnames(cor_matrix) <- assays
    p_matrix <- matrix(NA, nrow = n, ncol = n)
    rownames(p_matrix) <- assays
    colnames(p_matrix) <- assays
    for (i in 1:n) {
        for (j in 1:n) {
            if (i == j) {
                cor_matrix[i, j] <- 1
                p_matrix[i, j] <- 0
            }
            else if (i < j) {
                cat("Calculating:", assays[i], "vs", assays[j], "\n")
                cor_res <- multimodalallostery_calculate_correlation_simple(input_x = input_files[[assays[i]]], input_y = input_files[[assays[j]]], 
                  assay_x = assays[i], assay_y = assays[j], binding_sites_map = binding_sites_map)
                cor_matrix[i, j] <- cor_res$r
                cor_matrix[j, i] <- cor_res$r
                p_matrix[i, j] <- cor_res$p
                p_matrix[j, i] <- cor_res$p
            }
        }
    }
    return(list(cor_matrix = cor_matrix, p_matrix = p_matrix))
}

