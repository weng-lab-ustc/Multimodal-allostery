#' calculate correlation matrix
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_calculate_correlation_matrix()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param input_files Input or option used by this workflow; see Usage for the exact interface.
#' @param binder_order Input or option used by this workflow; see Usage for the exact interface.
#' @param binding_sites_map Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_calculate_correlation_matrix <- function (input_files, binder_order, binding_sites_map = list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), 
    RALGDS = c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67), PI3KCG = c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 
        63, 64, 70, 73), SOS1 = c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 
        66, 67, 70, 149, 153), K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74), K27 = c(21, 
        24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71), K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 
        97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138), K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 
        101, 102, 105, 107, 108, 125, 129, 133, 136, 137))) 
{
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

