#' calc or original
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_calc_or_original()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param df Input or option used by this workflow; see Usage for the exact interface.
#' @param region_residues Input or option used by this workflow; see Usage for the exact interface.
#' @param cat Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_calc_or_original <- function (df, region_residues, cat) 
{
    df[, `:=`(in_region, Pos_real %in% region_residues)]
    df[, `:=`(is_cat, category == cat)]
    a <- sum(df$in_region & df$is_cat)
    b <- sum(!df$in_region & df$is_cat)
    c <- sum(df$in_region & !df$is_cat)
    d <- sum(!df$in_region & !df$is_cat)
    tab <- matrix(c(a, b, c, d), nrow = 2)
    ft <- stats::fisher.test(tab)
    list(OR = unname(ft$estimate), p = ft$p.value, OR_low = ft$conf.int[1], OR_high = ft$conf.int[2])
}

