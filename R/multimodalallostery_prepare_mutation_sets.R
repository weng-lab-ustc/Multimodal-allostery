#' Prepare Mutation Sets.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param named_tables Value supplied for `named_tables`.
#' @param mutation_column Value supplied for `mutation_column`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_mutation_sets <- function(named_tables, mutation_column = "mutation") {
    lapply(named_tables, function(x) x[[mutation_column]])
}

