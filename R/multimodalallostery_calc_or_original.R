#' Calc or Original.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param df Input data frame or data table.
#' @param region_residues Value supplied for `region_residues`.
#' @param cat Value supplied for `cat`.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_calc_or_original <- function(df, region_residues, cat) {
    df[, `:=`(in_region, Pos_real %in% region_residues)]
    df[, `:=`(is_cat, category == cat)]
    a <- sum(df$in_region & df$is_cat)
    b <- sum(!df$in_region & df$is_cat)
    c <- sum(df$in_region & !df$is_cat)
    d <- sum(!df$in_region & !df$is_cat)
    tab <- matrix(c(a, b, c, d), nrow = 2)
    ft <- fisher.test(tab)
    list(OR = unname(ft$estimate), p = ft$p.value, OR_low = ft$conf.int[1], OR_high = ft$conf.int[2])
}

