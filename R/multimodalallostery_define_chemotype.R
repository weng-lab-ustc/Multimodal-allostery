#' Define Chemotype.
#'
#' Reusable classification/annotation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param aa Value supplied for `aa`.
#'
#' @return The input data with classifications or annotations added.
#' @export
multimodalallostery_define_chemotype <- function(aa) {
    aromatic <- c("F", "W", "Y")
    aliphatic <- c("A", "V", "I", "L", "M")
    polar_uncharged <- c("S", "T", "N", "Q", "C")
    positive <- c("K", "R", "H")
    negative <- c("D", "E")
    special <- c("G", "P")
    if (aa %in% aromatic) 
        return("Aromatic")
    if (aa %in% aliphatic) 
        return("Aliphatic")
    if (aa %in% polar_uncharged) 
        return("Polar uncharged")
    if (aa %in% positive) 
        return("Positive")
    if (aa %in% negative) 
        return("Negative")
    if (aa %in% special) 
        return("Special")
    return(NA)
}

