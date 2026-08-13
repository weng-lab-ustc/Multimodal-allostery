#' Plot Correlation Heatmap.
#'
#' Reusable plotting logic consolidated from the selected Figure/Panel implementations.
#'
#' @param cor_matrix Value supplied for `cor_matrix`.
#' @param p_matrix Value supplied for `p_matrix`.
#' @param output_file Optional output file path.
#' @param width Value supplied for `width`.
#' @param height Value supplied for `height`.
#'
#' @return A `ggplot2` plot object.
#' @export
multimodalallostery_plot_correlation_heatmap <- function(cor_matrix, p_matrix, output_file = NULL, width = 10, height = 8) {
    BI1_active <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55")
    BI1_inactive <- c("K27")
    BI2 <- c("K13", "K19")
    annotation_col <- data.frame(Group = ifelse(rownames(cor_matrix) %in% BI1_active, "BI1_active", ifelse(rownames(cor_matrix) %in% 
        BI2, "BI2", ifelse(rownames(cor_matrix) %in% BI1_inactive, "BI1_inactive", "Other"))))
    rownames(annotation_col) <- rownames(cor_matrix)
    annotation_colors <- list(Group = c(BI1_active = "#F4270C", BI1_inactive = "pink", BI2 = "#1B38A6"))
    heatmap_colors <- colorRampPalette(c("white", "#F4270C"))(100)
    sig_stars <- matrix("", nrow = nrow(p_matrix), ncol = ncol(p_matrix))
    for (i in 1:nrow(p_matrix)) {
        for (j in 1:ncol(p_matrix)) {
            if (i != j) {
                if (p_matrix[i, j] < 0.001) {
                  sig_stars[i, j] <- "***"
                }
                else if (p_matrix[i, j] < 0.01) {
                  sig_stars[i, j] <- "**"
                }
                else if (p_matrix[i, j] < 0.05) {
                  sig_stars[i, j] <- "*"
                }
                else {
                  sig_stars[i, j] <- "NS"
                }
            }
        }
    }
    diag(sig_stars) <- ""
    pheatmap(mat = cor_matrix, color = heatmap_colors, cluster_rows = FALSE, cluster_cols = FALSE, border_color = "grey90", 
        display_numbers = sig_stars, number_color = "black", number_size = 8, fontsize_number = 8, main = "Pearson Correlation of Mutational Effects (without BI sites)\nPer mutation", 
        xlab = "", ylab = "", annotation_col = annotation_col, annotation_row = annotation_col, annotation_colors = annotation_colors, 
        show_colnames = TRUE, show_rownames = TRUE, fontsize_row = 10, fontsize_col = 10, cellwidth = 40, cellheight = 40, 
        filename = output_file, width = width, height = height)
}

