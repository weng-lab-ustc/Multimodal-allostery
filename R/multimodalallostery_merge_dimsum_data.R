#' Merge Dimsum Data.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param merge_1 Value supplied for `merge_1`.
#' @param merge_2 Value supplied for `merge_2`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_merge_dimsum_data <- function(merge_1, merge_2) {
    a1 <- as.character(substitute(merge_1))
    a2 <- as.character(substitute(merge_2))
    merge_1[, `:=`(assay, a1)]
    merge_2[, `:=`(assay, a2)]
    output <- rbind(merge_1, merge_2)
    return(output)
}

