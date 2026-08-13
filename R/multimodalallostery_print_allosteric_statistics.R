#' Print Allosteric Statistics.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param data_plot Value supplied for `data_plot`.
#' @param allosteric_list Value supplied for `allosteric_list`.
#' @param assays Value supplied for `assays`.
#' @param reg_threshold Value supplied for `reg_threshold`.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_print_allosteric_statistics <- function(data_plot, allosteric_list, assays, reg_threshold) {
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("Statistical Report on Allosteric Sites\n")
    cat(rep("=", 80), "\n\n", sep = "")
    cat("\u3010Overall Statistics\u3011\n")
    cat(sprintf("threshold (reg_threshold): %.3f kcal/mol\n\n", reg_threshold))
    for (assayi in assays) {
        cat(sprintf("\u3010Assay: %s\u3011\n", assayi))
        cat(rep("-", 40), "\n", sep = "")
        data_current <- data_plot[assay == assayi, ]
        binding_sites <- allosteric_list[[assayi]][["Binding interface site"]]
        allosteric_gtp_sites <- allosteric_list[[assayi]][["Allosteric GTP pocket site"]]
        other_gtp_sites <- allosteric_list[[assayi]][["Other GTP pocket site"]]
        major_allosteric_sites <- allosteric_list[[assayi]][["Major allosteric site"]]
        cat("\n\u4F4D\u70B9\u6570\u91CF\u7EDF\u8BA1\uFF1A\n")
        cat(sprintf("  \u2022 Binding interface sites (\u7ED3\u5408\u754C\u9762\u4F4D\u70B9): %d\n", length(binding_sites)))
        cat(sprintf("  \u2022 Allosteric GTP pocket sites (\u53D8\u6784GTP\u53E3\u888B\u4F4D\u70B9): %d\n", length(allosteric_gtp_sites)))
        cat(sprintf("  \u2022 Other GTP pocket sites (\u5176\u4ED6GTP\u53E3\u888B\u4F4D\u70B9): %d\n", length(other_gtp_sites)))
        cat(sprintf("  \u2022 Major allosteric sites (\u4E3B\u8981\u53D8\u6784\u4F4D\u70B9): %d\n", length(major_allosteric_sites)))
        if (length(allosteric_gtp_sites) > 0) {
            cat("\nSpecific residues of the allosteric GTP-binding pocket\uFF1A\n")
            cat(sprintf("  Residue position: %s\n", paste(sort(allosteric_gtp_sites), collapse = ", ")))
            cat("\n  Detailed Information (mean ddG, distance, count):\n")
            for (pos in sort(allosteric_gtp_sites)) {
                pos_data <- data_current[Pos_real == pos, ]
                if (nrow(pos_data) > 0) {
                  cat(sprintf("    Pos %d: mean = %.3f, distance = %.2f \u00C5, count = %.1f\n", pos, pos_data$mean[1], pos_data$distance_bp[1], 
                    pos_data$count[1]))
                }
            }
        }
        if (length(major_allosteric_sites) > 0) {
            cat("\nSpecific residues at the primary allosteric site\uFF1A\n")
            cat(sprintf("  Residue position: %s\n", paste(sort(major_allosteric_sites), collapse = ", ")))
            cat("\n  Detailed Information (mean ddG, distance, count):\n")
            for (pos in sort(major_allosteric_sites)) {
                pos_data <- data_current[Pos_real == pos, ]
                if (nrow(pos_data) > 0) {
                  cat(sprintf("    Pos %d: mean = %.3f, distance = %.2f \u00C5, count = %.1f\n", pos, pos_data$mean[1], pos_data$distance_bp[1], 
                    pos_data$count[1]))
                }
            }
        }
        if (length(binding_sites) > 0) {
            cat("\nCombined with specific residues at the interface site\uFF1A\n")
            cat(sprintf("  Residue position: %s\n", paste(sort(binding_sites), collapse = ", ")))
        }
        if (length(other_gtp_sites) > 0) {
            cat("\nOther GTP pocket sites and specific residues\uFF1A\n")
            cat(sprintf("  Residue position: %s\n", paste(sort(other_gtp_sites), collapse = ", ")))
        }
        cat("\n")
    }
    cat(rep("=", 80), "\n", sep = "")
    cat("\u3010\u8DE8Assay\u6BD4\u8F83\u3011\n")
    cat(rep("-", 40), "\n", sep = "")
    all_allosteric_sites <- list()
    for (assayi in assays) {
        all_allosteric_sites[[assayi]] <- allosteric_list[[assayi]][["Allosteric GTP pocket site"]]
    }
    common_sites <- Reduce(intersect, all_allosteric_sites)
    if (length(common_sites) > 0) {
        cat("\nAllosteric GTP pocket site shared by all assays\uFF1A\n")
        cat(sprintf("  Residue position: %s\n", paste(sort(common_sites), collapse = ", ")))
    }
    else {
        cat("\nCommon allosteric GTP pocket sites across all assignments: None\n")
    }
    all_sites <- unique(unlist(all_allosteric_sites))
    cat(sprintf("\nTotal number of allosteric GTP pocket sites identified across all assays.: %d\n", length(all_sites)))
    if (length(all_sites) > 0) {
        cat(sprintf("  Full list: %s\n", paste(sort(all_sites), collapse = ", ")))
    }
    for (i in seq_along(assays)) {
        assayi <- assays[i]
        other_assays <- assays[-i]
        unique_sites <- setdiff(all_allosteric_sites[[assayi]], unlist(all_allosteric_sites[other_assays]))
        if (length(unique_sites) > 0) {
            cat(sprintf("\n%sUnique allosteric GTP pocket site\uFF1A\n", assayi))
            cat(sprintf("  Residue position: %s\n", paste(sort(unique_sites), collapse = ", ")))
        }
        else {
            cat(sprintf("\n%sUnique allosteric GTP pocket sites: None\n", assayi))
        }
    }
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("Statistical report concluded.\n")
    cat(rep("=", 80), "\n\n", sep = "")
}

