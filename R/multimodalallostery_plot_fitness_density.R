#' plot fitness density
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_plot_fitness_density()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param assay_type Input or option used by this workflow; see Usage for the exact interface.
#' @param block1 Input or option used by this workflow; see Usage for the exact interface.
#' @param block2 Input or option used by this workflow; see Usage for the exact interface.
#' @param block3 Input or option used by this workflow; see Usage for the exact interface.
#' @param output_file Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_plot_fitness_density <- function (assay_type, block1, block2, block3, output_file = NULL) 
{
    nor_fit <- multimodalallostery_normalize_fitness(block1 = block1, block2 = block2, block3 = block3)
    nor_fit_classified <- nor_fit %>% dplyr::mutate(mut_type = dplyr::case_when(Nham_aa == 0 & Nham_nt > 0 ~ "Synonymous", 
        STOP == TRUE | STOP_readthrough == TRUE ~ "Stop", Nham_aa > 0 & indel == FALSE & STOP == FALSE & STOP_readthrough == 
            FALSE ~ "Missense")) %>% dplyr::filter(!is.na(mut_type))
    cat("\nMutation distribution:", assay_type, "\n")
    print(table(nor_fit_classified$mut_type))
    nor_fit_plot <- nor_fit_classified %>% dplyr::mutate(mut_type = factor(mut_type, levels = c("Synonymous", "Missense", 
        "Stop")))
    p <- ggplot2::ggplot(nor_fit_plot, ggplot2::aes(x = fitness_normalized, color = mut_type)) + ggplot2::geom_density(linewidth = 1) + 
        ggplot2::scale_color_manual(values = c(Synonymous = "#1B38A6", Missense = "#F4AD0C", Stop = "#FF6A56")) + ggplot2::labs(title = paste0(toupper(assay_type), 
        " - Fitness Distribution"), x = "Normalized Fitness", y = "Density", color = "Mutation Type") + ggplot2::xlim(-1.5, 
        0.5) + ggplot2::theme_classic() + ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5, 
        size = 12), text = ggplot2::element_text(size = 10), legend.title = ggplot2::element_text(size = 10), legend.text = ggplot2::element_text(size = 9))
    if (!is.null(output_file)) {
        .ma_save_plot(output_file, p, width = 6, height = 4, units = "in")
        cat("Saved:", output_file, "\n")
    }
    return(p)
}

