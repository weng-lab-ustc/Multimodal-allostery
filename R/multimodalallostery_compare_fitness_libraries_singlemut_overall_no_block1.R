#' Compare Fitness Libraries Singlemut Overall No Block1.
#'
#' Reusable analysis/statistics logic consolidated from the selected Figure/Panel implementations.
#'
#' @param lib1_block2 Value supplied for `lib1_block2`.
#' @param lib1_block3 Value supplied for `lib1_block3`.
#' @param lib2_block2 Value supplied for `lib2_block2`.
#' @param lib2_block3 Value supplied for `lib2_block3`.
#' @param wt_aa Value supplied for `wt_aa`.
#' @param output_file Optional output file path.
#' @param x_lab Value supplied for `x_lab`.
#' @param y_lab Value supplied for `y_lab`.
#' @param main_title Value supplied for `main_title`.
#' @param point_alpha Value supplied for `point_alpha`.
#' @param plot_width Value supplied for `plot_width`.
#' @param plot_height Value supplied for `plot_height`.
#'
#' @return An analysis result, summary table, or numeric statistic.
#' @export
multimodalallostery_compare_fitness_libraries_singlemut_overall_no_block1 <- function(lib1_block2, lib1_block3, lib2_block2, lib2_block3, wt_aa, 
    output_file = NULL, x_lab = "Library 1 fitness", y_lab = "Library 2 fitness", main_title = "Comparison of fitness data between two libraries", 
    point_alpha = 0.3, plot_width = 5, plot_height = 4) {
    multimodalallostery_process_library_data <- function(block2, block3, wt_aa, suffix) {
        nor_fit <- nor_fitness(block2 = block2, block3 = block3)
        nor_fit_single <- nor_fitness_single_mut(input = nor_fit)
        nor_fit_single <- pos_id(nor_fit_single, wt_aa)
        fitness_data <- nor_fit_single[, c(1, 40, 41, 46, 48, 50, 52)]
        colnames(fitness_data) <- c("block", paste0("fitness", suffix), paste0("fitness_sigma", suffix), "Pos", "wtcodon", 
            "codon", "mt")
        return(fitness_data)
    }
    cat("Processing library 1 data...\n")
    fitness_data_1 <- multimodalallostery_process_library_data(lib1_block2, lib1_block3, wt_aa, "1")
    cat("Processing library 2 data...\n")
    fitness_data_2 <- multimodalallostery_process_library_data(lib2_block2, lib2_block3, wt_aa, "2")
    data <- merge(fitness_data_1, fitness_data_2, by = c("block", "Pos", "wtcodon", "codon", "mt"), all = FALSE)
    setDT(data)
    cat(paste("Total variants after merging:", nrow(data), "\n"))
    complete_cases <- complete.cases(data$fitness1, data$fitness2)
    data_complete <- data[complete_cases, ]
    if (nrow(data_complete) < 2) {
        warning("Insufficient complete cases for plot")
        return(ggplot() + labs(title = main_title, x = x_lab, y = y_lab) + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") + 
            theme_minimal())
    }
    cor_test <- cor.test(data_complete$fitness1, data_complete$fitness2, method = "pearson", use = "complete.obs")
    r_value <- round(cor_test$estimate, 3)
    p_value <- round(cor_test$p.value, 4)
    p <- ggplot(data_complete, aes(x = fitness1, y = fitness2)) + geom_point(color = "#75C2F6", alpha = point_alpha, size = 1.5) + 
        coord_cartesian(xlim = c(-1.5, 1), ylim = c(-1.5, 0.7)) + labs(title = main_title, x = x_lab, y = y_lab) + annotate("text", 
        x = -1.4, y = 0.45, label = paste0("r = ", r_value, "\np = ", ifelse(p_value < 1e-04, "< 0.0001", p_value)), hjust = 0, 
        vjust = 1, size = 3, color = "black") + theme_classic(base_size = 10) + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, 
        size = 11), axis.text = element_text(size = 10, colour = "black"), axis.text.x = element_text(angle = 90, vjust = 0.5, 
        hjust = 1), axis.title = element_text(size = 10), legend.position = "none", plot.margin = margin(10, 10, 10, 10), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
    print(p)
    if (!is.null(output_file)) {
        ggsave(filename = output_file, plot = p, device = cairo_pdf, width = plot_width, height = plot_height, units = "in", 
            dpi = 300)
        cat(paste("Plot saved to:", output_file, "\n"))
    }
    return(list(data = data, plot = p))
}

