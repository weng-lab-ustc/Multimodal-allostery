#' Print Mutation List.
#'
#' Reusable utility logic consolidated from the selected Figure/Panel implementations.
#'
#' @param mutations Value supplied for `mutations`.
#' @param title Value supplied for `title`.
#'
#' @return The transformed result produced by the function.
#' @export
multimodalallostery_print_mutation_list <- function(mutations, title) {
    cat("\n", title, "\n", sep = "")
    cat("\u7A81\u53D8\u6570:", length(mutations), "\n")
    if (length(mutations) > 0) {
        for (i in seq(1, length(mutations), by = 10)) {
            end_idx <- min(i + 9, length(mutations))
            cat(paste(mutations[i:end_idx], collapse = ", "), "\n")
        }
    }
    else {
        cat("None\n")
    }
}

