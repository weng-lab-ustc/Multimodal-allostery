#' Prepare Allosteric Pair Counts.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data Input data frame or data table.
#' @param pair_order Value supplied for `pair_order`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_allosteric_pair_counts <- function(data, pair_order) {
    dplyr::mutate(tidyr::pivot_longer(dplyr::select(dplyr::mutate(dplyr::filter(data, Pair %in% pair_order), Pair = factor(Pair, 
        levels = pair_order)), Pair, Allosteric_only_in_X, Allosteric_only_in_Y), cols = c(Allosteric_only_in_X, Allosteric_only_in_Y), 
        names_to = "Type", values_to = "Count"), Type = factor(Type, levels = c("Allosteric_only_in_X", "Allosteric_only_in_Y"), 
        labels = c("Allosteric only in X", "Allosteric only in Y")))
}

