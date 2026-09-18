#' print allosteric statistics
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_print_allosteric_statistics()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param data_plot Input or option used by this workflow; see Usage for the exact interface.
#' @param allosteric_list Input or option used by this workflow; see Usage for the exact interface.
#' @param assays Input or option used by this workflow; see Usage for the exact interface.
#' @param reg_threshold Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_print_allosteric_statistics <- function (data_plot, allosteric_list, assays, reg_threshold) 
{
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("\u53d8\u6784\u4f4d\u70b9\u7edf\u8ba1\u62a5\u544a\n")
    cat(rep("=", 80), "\n\n", sep = "")
    cat("\u3010\u603b\u4f53\u7edf\u8ba1\u3011\n")
    cat(sprintf("\u56de\u5f52\u9608\u503c (reg_threshold): %.3f kcal/mol\n", reg_threshold))
    cat(sprintf("\u57fa\u4e8e\u6240\u6709 %d \u4e2abinder\u7684binding sites\u8ba1\u7b97\n\n", length(assays)))
    for (assayi in assays) {
        cat(sprintf("\u3010Assay: %s\u3011\n", assayi))
        cat(rep("-", 40), "\n", sep = "")
        data_current <- data_plot[assay == assayi, ]
        binding_sites <- allosteric_list[[assayi]][["Binding interface site"]]
        allosteric_gtp_sites <- allosteric_list[[assayi]][["Allosteric GTP pocket site"]]
        other_gtp_sites <- allosteric_list[[assayi]][["Other GTP pocket site"]]
        major_allosteric_sites <- allosteric_list[[assayi]][["Major allosteric site"]]
        cat("\n\u4f4d\u70b9\u6570\u91cf\u7edf\u8ba1\uff1a\n")
        cat(sprintf("  \u2022 Binding interface sites (\u7ed3\u5408\u754c\u9762\u4f4d\u70b9): %d\n", length(binding_sites)))
        cat(sprintf("  \u2022 Allosteric GTP pocket sites (\u53d8\u6784GTP\u53e3\u888b\u4f4d\u70b9): %d\n", length(allosteric_gtp_sites)))
        cat(sprintf("  \u2022 Other GTP pocket sites (\u5176\u4ed6GTP\u53e3\u888b\u4f4d\u70b9): %d\n", length(other_gtp_sites)))
        cat(sprintf("  \u2022 Major allosteric sites (\u4e3b\u8981\u53d8\u6784\u4f4d\u70b9): %d\n", length(major_allosteric_sites)))
        if (length(allosteric_gtp_sites) > 0) {
            cat("\n\u53d8\u6784GTP\u53e3\u888b\u4f4d\u70b9\u5177\u4f53\u6b8b\u57fa\uff1a\n")
            cat(sprintf("  \u6b8b\u57fa\u4f4d\u7f6e: %s\n", paste(sort(allosteric_gtp_sites), collapse = ", ")))
            cat("\n  \u8be6\u7ec6\u4fe1\u606f (mean ddG, distance, count):\n")
            for (pos in sort(allosteric_gtp_sites)) {
                pos_data <- data_current[Pos_real == pos, ]
                if (nrow(pos_data) > 0) {
                  cat(sprintf("    Pos %d: mean = %.3f, distance = %.2f \u00c5, count = %.1f\n", pos, pos_data$mean[1], pos_data$distance_bp[1], 
                    pos_data$count[1]))
                }
            }
        }
        if (length(major_allosteric_sites) > 0) {
            cat("\n\u4e3b\u8981\u53d8\u6784\u4f4d\u70b9\u5177\u4f53\u6b8b\u57fa\uff1a\n")
            cat(sprintf("  \u6b8b\u57fa\u4f4d\u7f6e: %s\n", paste(sort(major_allosteric_sites), collapse = ", ")))
            cat("\n  \u8be6\u7ec6\u4fe1\u606f (mean ddG, distance, count):\n")
            for (pos in sort(major_allosteric_sites)) {
                pos_data <- data_current[Pos_real == pos, ]
                if (nrow(pos_data) > 0) {
                  cat(sprintf("    Pos %d: mean = %.3f, distance = %.2f \u00c5, count = %.1f\n", pos, pos_data$mean[1], pos_data$distance_bp[1], 
                    pos_data$count[1]))
                }
            }
        }
        if (length(binding_sites) > 0) {
            cat("\n\u7ed3\u5408\u754c\u9762\u4f4d\u70b9\u5177\u4f53\u6b8b\u57fa\uff1a\n")
            cat(sprintf("  \u6b8b\u57fa\u4f4d\u7f6e: %s\n", paste(sort(binding_sites), collapse = ", ")))
        }
        if (length(other_gtp_sites) > 0) {
            cat("\n\u5176\u4ed6GTP\u53e3\u888b\u4f4d\u70b9\u5177\u4f53\u6b8b\u57fa\uff1a\n")
            cat(sprintf("  \u6b8b\u57fa\u4f4d\u7f6e: %s\n", paste(sort(other_gtp_sites), collapse = ", ")))
        }
        cat("\n")
    }
    cat(rep("=", 80), "\n", sep = "")
    cat("\u3010\u8de8Assay\u6bd4\u8f83\u3011\n")
    cat(rep("-", 40), "\n", sep = "")
    all_allosteric_sites <- list()
    for (assayi in assays) {
        all_allosteric_sites[[assayi]] <- allosteric_list[[assayi]][["Allosteric GTP pocket site"]]
    }
    common_sites <- Reduce(intersect, all_allosteric_sites)
    if (length(common_sites) > 0) {
        cat("\n\u6240\u6709assay\u5171\u6709\u7684\u53d8\u6784GTP\u53e3\u888b\u4f4d\u70b9\uff1a\n")
        cat(sprintf("  \u6b8b\u57fa\u4f4d\u7f6e: %s\n", paste(sort(common_sites), collapse = ", ")))
    }
    else {
        cat("\n\u6240\u6709assay\u5171\u6709\u7684\u53d8\u6784GTP\u53e3\u888b\u4f4d\u70b9\uff1a\u65e0\n")
    }
    all_sites <- unique(unlist(all_allosteric_sites))
    cat(sprintf("\n\u6240\u6709assay\u4e2d\u51fa\u73b0\u7684\u53d8\u6784GTP\u53e3\u888b\u4f4d\u70b9\u603b\u6570: %d\n", length(all_sites)))
    if (length(all_sites) > 0) {
        cat(sprintf("  \u5b8c\u6574\u5217\u8868: %s\n", paste(sort(all_sites), collapse = ", ")))
    }
    for (i in seq_along(assays)) {
        assayi <- assays[i]
        other_assays <- assays[-i]
        unique_sites <- setdiff(all_allosteric_sites[[assayi]], unlist(all_allosteric_sites[other_assays]))
        if (length(unique_sites) > 0) {
            cat(sprintf("\n%s\u7279\u6709\u7684\u53d8\u6784GTP\u53e3\u888b\u4f4d\u70b9\uff1a\n", assayi))
            cat(sprintf("  \u6b8b\u57fa\u4f4d\u7f6e: %s\n", paste(sort(unique_sites), collapse = ", ")))
        }
        else {
            cat(sprintf("\n%s\u7279\u6709\u7684\u53d8\u6784GTP\u53e3\u888b\u4f4d\u70b9\uff1a\u65e0\n", assayi))
        }
    }
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("\u7edf\u8ba1\u62a5\u544a\u7ed3\u675f\n")
    cat(rep("=", 80), "\n\n", sep = "")
}

