#' Calculate Threshold.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data Input data frame or data table.
#' @param assay_sele Selected assay identifier.
#' @param anno Annotation data or annotation file path used by the analysis.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_calculate_threshold <- function(data, assay_sele, anno) {
    anno_data <- fread(anno)
    site_ddG <- data[, .(mean_ddG = sum(abs(ddG)/ddG_std^2, na.rm = TRUE)/sum(1/ddG_std^2, na.rm = TRUE)), by = Pos_real]
    site_sigma <- data[, .(sigma = sqrt(1/sum(1/ddG_std^2, na.rm = TRUE))), by = Pos_real]
    site_stats <- merge(site_ddG, site_sigma, by = "Pos_real")
    site_stats <- merge(site_stats, anno_data, by.x = "Pos_real", by.y = "Pos", all.x = TRUE)
    scHAmin_col <- paste0("scHAmin_ligand_", assay_sele)
    binding_sites <- site_stats[get(scHAmin_col) < 5, Pos_real]
    threshold <- site_stats[Pos_real %in% binding_sites, sum(abs(mean_ddG)/sigma^2, na.rm = TRUE)/sum(1/sigma^2, na.rm = TRUE)]
    return(threshold)
}

