library(krasddpcams)
library(data.table)
library(wlab.block)
#options(timeout = 600)
#remotes::install_github(
#  "weng-lab-ustc/Multimodal-allostery",
#  dependencies = TRUE
#)
library(multimodalallostery)


setwd("~/Downloads/Multimodal-allostery-main/Supplementary_data/")

# ✅
# ---- Figure1E ----
run_figure1_e <- function() {
  result <- local({
    Abundance <- multimodalallostery_normalize_fitness(block1 = "./fitness_RData/CW_RAS_abundance_1_fitness_replicates_fullseq.RData", 
                                   block2 = "./fitness_RData_merge_version/MA_RAS_abundance_2_fitness_replicates_fullseq_merge.RData", 
                                   block3 = "./fitness_RData_merge_version/MA_RAS_abundance_3_fitness_replicates_fullseq_merge.RData")
    Abundance <- multimodalallostery_classify_mutation_type(Abundance)
    cat("\nMutation distribution: Abundance\n")
    print(table(Abundance$mut_type))
    Abundance_plot <- multimodalallostery_plot_fitness_density(Abundance, assay_type = "Abundance")
    ggplot2::ggsave("./results/f1e_Abundance_normalized_density.pdf", 
           Abundance_plot, width = 6, height = 4, units = "in")
    K13 <- multimodalallostery_normalize_fitness(block1 = "./fitness_RData_merge_version/MA_RAS_binding_K13_1_fitness_replicates_fullseq_merge.RData", 
                             block2 = "./fitness_RData_merge_version/MA_RAS_binding_K13_2_fitness_replicates_fullseq_merge.RData", 
                             block3 = "./fitness_RData_merge_version/MA_RAS_binding_K13_3_fitness_replicates_fullseq_merge.RData")
    K13 <- multimodalallostery_classify_mutation_type(K13)
    cat("\nMutation distribution: K13\n")
    print(table(K13$mut_type))
    K13_plot <- multimodalallostery_plot_fitness_density(K13, assay_type = "K13")
    ggplot2::ggsave("./results/f1e_K13_normalized_density.pdf", 
           K13_plot, width = 6, height = 4, units = "in")
    K19 <- multimodalallostery_normalize_fitness(block1 = "./fitness_RData_merge_version/MA_RAS_binding_K19_1_fitness_replicates_fullseq_merge.RData", 
                             block2 = "./fitness_RData_merge_version/MA_RAS_binding_K19_2_fitness_replicates_fullseq_merge.RData", 
                             block3 = "./fitness_RData_merge_version/MA_RAS_binding_K19_3_fitness_replicates_fullseq_merge.RData")
    K19 <- multimodalallostery_classify_mutation_type(K19)
    cat("\nMutation distribution: K19\n")
    print(table(K19$mut_type))
    K19_plot <- multimodalallostery_plot_fitness_density(K19, assay_type = "K19")
    ggplot2::ggsave("./results/f1e_K19_normalized_density.pdf", 
           K19_plot, width = 6, height = 4, units = "in")
  })
  invisible(result)
}


run_figure1_e()

# ✅
# ---- Figure1F ----
run_figure1_f <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    stability_nor_df <- krasddpcams::krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./fitness_RData/CW_RAS_abundance_1_fitness_replicates_fullseq.RData", 
                                                                               block2_dimsum_df = "./fitness_RData/MA_RAS_abundance_2_fitness_replicates_fullseq.RData", 
                                                                               block3_dimsum_df = "./fitness_RData/MA_RAS_abundance_3_fitness_replicates_fullseq.RData")
    K13_nor_df <- krasddpcams::krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./fitness_RData/MA_RAS_binding_K13_1_fitness_replicates_fullseq.RData", 
                                                                         block2_dimsum_df = "./fitness_RData/MA_RAS_binding_K13_2_fitness_replicates_fullseq.RData", 
                                                                         block3_dimsum_df = "./fitness_RData/MA_RAS_binding_K13_3_fitness_replicates_fullseq.RData")
    K19_nor_df <- krasddpcams::krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./fitness_RData/MA_RAS_binding_K19_1_fitness_replicates_fullseq.RData", 
                                                                         block2_dimsum_df = "./fitness_RData/MA_RAS_binding_K19_2_fitness_replicates_fullseq.RData", 
                                                                         block3_dimsum_df = "./fitness_RData/MA_RAS_binding_K19_3_fitness_replicates_fullseq.RData")
    RAF1_nor_df <- krasddpcams::krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./fitness_RData/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData", 
                                                                          block2_dimsum_df = "./fitness_RData/MA_RAS_binding_RAF_2_fitness_replicates_fullseq.RData", 
                                                                          block3_dimsum_df = "./fitness_RData/MA_RAS_binding_RAF_3_fitness_replicates_fullseq.RData")
    stab <- stability_nor_df
    K13 <- K13_nor_df
    K19 <- K19_nor_df
    RAF1 <- RAF1_nor_df
    all_data_K13 <- multimodalallostery_merge_dimsum_data(stab, K13)
    all_data_K19 <- multimodalallostery_merge_dimsum_data(stab, K19)
    all_data_RAF1 <- multimodalallostery_merge_dimsum_data(stab, RAF1)
    all_data_pos_K13 <- multimodalallostery_identify_mutation_positions(all_data_K13, wt_aa)
    all_data_pos_K19 <- multimodalallostery_identify_mutation_positions(all_data_K19, wt_aa)
    all_data_pos_RAF1 <- multimodalallostery_identify_mutation_positions(all_data_RAF1, wt_aa)
    anno <- data.table::fread("./anno_final_for_5.csv")
    anno[, `:=`(Pos_real, Pos)]
    summery <- wlab.block::ddG_data_assay(input = "./weights_Binding_RAF.txt", 
                                          wt_aa = wt_aa)
    anno <- merge(anno, summery, by = "Pos_real", all = T)
    plot_K13 <- multimodalallostery_plot_binding_fitness(input = all_data_pos_K13, assay_sele = "K13", anno = anno)
    plot_K19 <- multimodalallostery_plot_binding_fitness(input = all_data_pos_K19, assay_sele = "K19", anno = anno)
    plot_RAF1 <- multimodalallostery_plot_binding_fitness(input = all_data_pos_RAF1, assay_sele = "RAF1", anno = anno)
    print(plot_K13)
    print(plot_K19)
    print(plot_RAF1)
    ggplot2::ggsave("./results/f1f_K13_binding_vs_abundance_fitness.pdf", 
                    plot = plot_K13, device = grDevices::cairo_pdf, height = 4, width = 4)
    ggplot2::ggsave("./results/f1f_K19_binding_vs_abundance_fitness.pdf", 
                    plot = plot_K19, device = grDevices::cairo_pdf, height = 4, width = 4)
    ggplot2::ggsave("./results/f1f_RAF1_binding_vs_abundance_fitness.pdf", 
                    plot = plot_RAF1, device = grDevices::cairo_pdf, height = 4, width = 4)
  })
  invisible(result)
}

run_figure1_f()



# ✅
# ---- Figure1G ----
run_figure1_g <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    anno <- data.table::fread("./anno_final_for_5.csv")
    p_k13 <- multimodalallostery_plot_scatter_dd_gb_dd_gf(ddG1 = "./weights_Folding.txt", 
                                      assay1 = "folding", ddG2 = "./weights_Binding_K13.txt", 
                                      assay2 = "K13", anno = anno, binder = "K13")
    p_k13
    ggplot2::ggsave("./results/f1g_scatter_ddGb_ddGf_K13.pdf", 
                    p_k13, device = grDevices::cairo_pdf, height = 4, width = 4)
    p_k19 <- multimodalallostery_plot_scatter_dd_gb_dd_gf(ddG1 = "./weights_Folding.txt", 
                                      assay1 = "folding", ddG2 = "./weights_Binding_K19.txt", 
                                      assay2 = "K19", anno = anno, binder = "K19")
    p_k19
    ggplot2::ggsave("./results/f1g_scatter_ddGb_ddGf_K19.pdf", 
                    p_k19, device = "pdf", height = 4, width = 4)
    p_RAF1 <- multimodalallostery_plot_scatter_dd_gb_dd_gf(ddG1 = "./weights_Folding.txt", 
                                       assay1 = "folding", ddG2 = "./weights_Binding_RAF.txt", 
                                       assay2 = "RAF1", anno = anno, binder = "RAF1")
    p_RAF1
    ggplot2::ggsave("./results/f1g_scatter_ddGb_ddGf_RAF1.pdf", 
                    p_RAF1, device = "pdf", height = 4, width = 4)
  })
  invisible(result)
}

run_figure1_g()



# ✅
# ---- Figure1H ----
run_figure1_h <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    multimodalallostery_dd_g_heatmap(input = "./weights_Folding.txt", 
                 wt_aa = wt_aa, title = "KRAS-Folding free energy changes", legend_limits = c(-1.5, 3.3))
    ggplot2::ggsave("./results/f1h_KRAS-Folding free energy changes heatmap.pdf", 
                    device = "pdf", height = 6, width = 20)
    multimodalallostery_dd_g_heatmap(input = "./weights_Binding_RAF.txt", 
                 wt_aa = wt_aa, title = "KRAS-RAF1 binding free energy changes", legend_limits = c(-1.5, 3.3))
    ggplot2::ggsave("./results/f1h_KRAS-RAF1 binding free energy changes.pdf", 
                    device = "pdf", height = 6, width = 20)
    multimodalallostery_dd_g_heatmap(input = "./weights_Binding_SOS.txt", 
                 wt_aa = wt_aa, title = "KRAS-SOS1 binding free energy changes heatmap", legend_limits = c(-1.5, 3.3))
    ggplot2::ggsave("./results/f1h_KRAS-SOS1 binding free energy changes.pdf", 
                    device = "pdf", height = 6, width = 20)
    multimodalallostery_dd_g_heatmap(input = "./weights_Binding_K55.txt", 
                 wt_aa = wt_aa, title = "KRAS-K55 binding free energy changes heatmap", legend_limits = c(-1.5, 3.3))
    ggplot2::ggsave("./results/f1h_KRAS-K55 binding free energy changes.pdf", 
                    device = "pdf", height = 6, width = 20)
    multimodalallostery_dd_g_heatmap(input = "./weights_Binding_K27.txt", 
                 wt_aa = wt_aa, title = "KRAS-K27 binding free energy changes heatmap", legend_limits = c(-1.5, 3.3))
    ggplot2::ggsave("./results/f1h_KRAS-K27 binding free energy changes.pdf", 
                    device = "pdf", height = 6, width = 20)
    multimodalallostery_dd_g_heatmap(input = "./weights_Binding_RAL.txt", 
                 wt_aa = wt_aa, title = "KRAS-RALGDS binding free energy changes heatmap", legend_limits = c(-1.5, 3.3))
    ggplot2::ggsave("./results/f1h_KRAS-RALGDS binding free energy changes.pdf", 
                    device = "pdf", height = 6, width = 20)
    multimodalallostery_dd_g_heatmap(input = "./weights_Binding_PI3.txt", 
                 wt_aa = wt_aa, title = "KRAS-PIK3CG binding free energy changes heatmap", legend_limits = c(-1.5, 3.3))
    ggplot2::ggsave("./results/f1h_KRAS-PIK3CG binding free energy changes.pdf", 
                    device = "pdf", height = 6, width = 20)
    multimodalallostery_dd_g_heatmap(input = "./weights_Binding_K13.txt", 
                 wt_aa = wt_aa, title = "KRAS-DARPin_K13 binding free energy changes heatmap", legend_limits = c(-1.5, 3.3))
    ggplot2::ggsave("./results/f1h_K13 heatmap.pdf", device = "pdf", 
                    height = 6, width = 20)
    multimodalallostery_dd_g_heatmap(input = "./weights_Binding_K19.txt", 
                 wt_aa = wt_aa, title = "KRAS-DARPin_K19 binding free energy changes", legend_limits = c(-1.5, 3.3))
    ggplot2::ggsave("./results/f1h_KRAS-DARPin_K19 binding free energy changes.pdf", 
                    device = "pdf", height = 6, width = 20)
  })
  invisible(result)
}

run_figure1_h()


# ✅
# ---- Figure1I ----
run_figure1_i <- function() {
  result <- local({
    KRAS_sequence <- "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    K13_binding_interface_site <- c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 
                                    137, 138)
    RAF1_binding_interface_site <- c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71)
    GTP_pocket <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    functional_loop <- data.frame(xstart = c(10, 25, 58), xend = c(17, 40, 76), col = c("P-loop", "switch I", "switch II"))
    common_core_residues <- c(4, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 40, 42, 44, 46, 51, 52, 53, 
                              54, 55, 56, 57, 58, 68, 71, 72, 75, 77, 78, 79, 80, 81, 82, 83, 84, 89, 90, 92, 93, 96, 97, 99, 100, 101, 103, 109, 110, 
                              111, 112, 113, 114, 115, 116, 118, 125, 130, 133, 134, 137, 139, 141, 142, 143, 144, 145, 146, 151, 152, 155, 156, 157, 
                              158, 159, 160, 162, 163)
    rects_sheet <- data.frame(xstart = c(3, 38, 51, 77, 109, 139), xend = c(9, 44, 57, 84, 115, 143), col = c("sheet1", "sheet2", 
                                                                                                              "sheet3", "sheet4", "sheet5", "sheet6"))
    rects_helix <- data.frame(xstart = c(15, 67, 87, 127, 148), xend = c(24, 73, 104, 136, 166), col = c("helix1", "helix2", 
                                                                                                         "helix3", "helix4", "helix5"))
    sequence_annotation <- multimodalallostery_prepare_sequence_annotation(sequence = KRAS_sequence, binding_sites_1 = K13_binding_interface_site, 
                                                       binding_sites_2 = RAF1_binding_interface_site)
    p <- multimodalallostery_plot_sequence_annotation(prepared_sequence = sequence_annotation, gtp_pocket = GTP_pocket, functional_loops = functional_loop, 
                                  core_residues = common_core_residues, beta_sheets = rects_sheet, alpha_helices = rects_helix)
    print(p)
    ggplot2::ggsave("./results/f1i_sequence_structure_rotated.pdf", 
                    plot = p, device = "pdf", height = 4, width = 20, dpi = 300)
  })
  invisible(result)
}


run_figure1_i()


# ✅
# ---- Figure2A ----
run_figure2_a <- function() {
  result <- local({
    ddG_file <- "./weights_Binding_K13.txt"
    binding_interface_residues <- c(63, 105, 106, 98, 107, 101, 102, 99, 136, 95, 137, 94, 133, 90, 129, 87, 91, 88)
    position_labels <- c(`106` = "S106", `105` = "D105", `63` = "E63", `98` = "E98", `107` = "E107", `101` = "K101", `102` = "R102", 
                         `99` = "Q99", `136` = "S136", `95` = "H95", `137` = "Y137", `94` = "H94", `133` = "L133", `90` = "F90", `129` = "Q129", 
                         `87` = "T87", `91` = "E91", `88` = "K88")
    p <- multimodalallostery_plot_binding_interface_residue_median_dd_g_heatmap(ddG_file = ddG_file, binding_sites = binding_interface_residues, 
                                                            position_labels = position_labels, title = "K13 Binding Interface Residues - Median ΔΔGb")
    print(p)
    ggplot2::ggsave("./results/f2a_K13_BI_residues_median_ddG_heatmap.pdf", 
                    p, device = "pdf", height = 8, width = 3)
    ddG_file <- "./weights_Binding_K19.txt"
    binding_interface_residues <- c(98, 107, 101, 102, 99, 136, 95, 137, 94, 133, 90, 129, 87, 91, 88, 68, 108)
    position_labels <- c(`98` = "E98", `107` = "E107", `101` = "K101", `102` = "R102", `99` = "Q99", `136` = "S136", `95` = "H95", 
                         `137` = "Y137", `94` = "H94", `133` = "L133", `90` = "F90", `129` = "Q129", `87` = "T87", `91` = "E91", `88` = "K88", 
                         `68` = "R68", `108` = "D108")
    p <- multimodalallostery_plot_binding_interface_residue_median_dd_g_heatmap(ddG_file = ddG_file, binding_sites = binding_interface_residues, 
                                                            position_labels = position_labels, title = "K19 Binding Interface Residues - Median ΔΔGb")
    print(p)
    ggplot2::ggsave("./results/f2a_K19_BI_residues_median_ddG_heatmap.pdf", 
                    p, device = "pdf", height = 7, width = 3)
  })
  invisible(result)
}



run_figure2_a()


# ✅
# ---- Figure2C ----
run_figure2_c <- function() {
  result <- local({
    residues_list <- c("K88", "E91", "T87", "Q129", "F90", "L133", "H94", "Y137", "H95", "R68", "S136", "Q99", "R102", "K101", 
                       "E107", "E98")
    multimodalallostery_plot_dd_g_beeswarm(ddG_file = "./weights_Binding_K13.txt", 
                       assay_sele = "K13", residues = residues_list, output_file = "./results/f2c_Figure_sites_mut_ddG_beeswarm_K13.pdf")
    multimodalallostery_plot_dd_g_beeswarm(ddG_file = "./weights_Binding_K19.txt", 
                       assay_sele = "K19", residues = residues_list, output_file = "./results/f2c_Figure_sites_mut_ddG_beeswarm_K19.pdf")
  })
  invisible(result)
}

run_figure2_c()



# ✅
# ---- Figure2D ----
run_figure2_d <- function() {
  result <- local({
    ddG1 <- krasddpcams::krasddpcams__read_ddG(ddG = "./weights_Binding_K13.txt", 
                                               assay_sele = "K13")
    ddG2 <- krasddpcams::krasddpcams__read_ddG(ddG = "./weights_Binding_K19.txt", 
                                               assay_sele = "K19")
    interface_sites <- c(98, 107, 101, 102, 99, 136, 95, 137, 94, 133, 90, 129, 87, 91, 88)
    comparison_data <- multimodalallostery_prepare_interface_ddg_comparison(data_x = ddG1, data_y = ddG2, interface_x = interface_sites, interface_y = interface_sites, 
                                                        assay_x = "K13", assay_y = "K19")
    analysis <- multimodalallostery_analyze_correlation_outliers(comparison_data, x_var = "K13", y_var = "K19", num_outliers = 10)
    cat("=== Correlation Statistics ===\n")
    cat("Pearson r =", round(analysis$correlation$estimate, 4), "\n")
    cat("p-value =", format(analysis$correlation$p.value, scientific = TRUE, digits = 4), "\n")
    cat("Sample size (n) =", analysis$correlation$parameter + 2, "\n")
    print(analysis$outliers)
    cat("Number of outliers that might be due to measurement noise:", analysis$noise_count, "/", analysis$num_outliers, "\n")
    p <- multimodalallostery_plot_ddg_correlation_outliers(analysis, x_label = "ΔΔG for K13 Binding Interface Mutations (kcal/mol)", y_label = "ΔΔG for K19 Binding Interface Mutations (kcal/mol)")
    print(p)
    ggplot2::ggsave("./results/f2d_scatter_plot_compare_K13_K19_binding_interface_label_outlier_per_mutations.pdf", 
                    plot = p, device = "pdf", height = 6, width = 6, dpi = 300)
  })
  invisible(result)
}


run_figure2_d()


# ✅
# ---- Figure3A ----
run_figure3_a <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    rects_sheet <- data.frame(xstart = c(3, 38, 51, 77, 109, 139), xend = c(9, 44, 57, 84, 115, 143), col = c("β1", "β2", "β3", 
                                                                                                              "β4", "β5", "β6"))
    rects_helix <- data.frame(xstart = c(15, 67, 87, 127, 148), xend = c(24, 73, 104, 136, 166), col = c("α1", "α2", "α3", 
                                                                                                         "α4", "α5"))
    multimodalallostery_manhatta_plot_single_assay(input = "./weights_Binding_K13.txt", 
                               assay_sele = "K13", anno = "./anno_final_for_5.csv", 
                               rects_sheet = rects_sheet, rects_alpha = rects_helix, wt_aa = wt_aa)
    ggplot2::ggsave("./results/f3a_K13_Manhattan_plot.pdf", 
                    device = grDevices::cairo_pdf, height = 9, width = 12)
    multimodalallostery_manhatta_plot_single_assay(input = "./weights_Binding_RAF.txt", 
                               assay_sele = "RAF1", anno = "./anno_final_for_5.csv", 
                               rects_sheet = rects_sheet, rects_alpha = rects_helix, wt_aa = wt_aa)
    ggplot2::ggsave("./results/f3a_RAF1_Manhattan_plot.pdf", 
                    device = grDevices::cairo_pdf, height = 9, width = 12)
  })
  invisible(result)
}


run_figure3_a()




# ✅
# ---- Figure3C ----
run_figure3_c <- function() {
  result <- local({
    plot_raf1 <- multimodalallostery_plot_energy_distance_decay_expfit_annotation(input = "./weights_Binding_RAF.txt", 
                                                              assay_sele = "RAF1", anno_file = "./anno_final_for_5.csv")
    print(plot_raf1)
    ggplot2::ggsave("./results/f3c_RAF1_distance_decay all mutations.pdf", 
                    plot_raf1, device = grDevices::cairo_pdf, width = 2.5, height = 2.5, units = "in", dpi = 300)
    plot_k13 <- multimodalallostery_plot_energy_distance_decay_expfit_annotation(input = "./weights_Binding_K13.txt", 
                                                             assay_sele = "K13", anno_file = "./anno_final_for_5.csv")
    print(plot_k13)
    ggplot2::ggsave("./results/f3c_K13_distance_decay all mutations.pdf", 
                    plot_k13, device = grDevices::cairo_pdf, width = 2.5, height = 2.5, units = "in", dpi = 300)
  })
  invisible(result)
}

run_figure3_c()



# ✅
# ---- Figure3D ----
run_figure3_d <- function() {
  result <- local({
    plot_raf1 <- multimodalallostery_plot_energy_distance_decay_directional_no_filter(input = "./weights_Binding_RAF.txt", 
                                                                  assay_sele = "RAF1", anno_file = "./anno_final_for_5.csv")
    print(plot_raf1)
    ggplot2::ggsave("./results/f3d_RAF1_directional_decay.pdf", 
                    plot_raf1, device = grDevices::cairo_pdf, width = 2.5, height = 2.5, units = "in", dpi = 300)
    plot_K13 <- multimodalallostery_plot_energy_distance_decay_directional_no_filter(input = "./weights_Binding_K13.txt", 
                                                                 assay_sele = "K13", anno_file = "./anno_final_for_5.csv")
    print(plot_K13)
    ggplot2::ggsave("./results/f3d_K13_directional_decay.pdf", 
                    plot_K13, device = grDevices::cairo_pdf, width = 2.5, height = 2.5, units = "in", dpi = 300)
  })
  invisible(result)
}

run_figure3_d()


# ✅
# ---- Figure3E ----
run_figure3_e <- function() {
  result <- local({
    contact_shell <- data.table::fread("./5binder_contact_shell2.csv")
    plot_raf1 <- multimodalallostery_plot_energy_distance_decay_expfit_contact_shell(input = "./weights_Binding_RAF.txt", 
                                                                 assay_sele = "RAF1", contact_shell = contact_shell)
    print(plot_raf1)
    ggplot2::ggsave("./results/f3e_RAF1_contact_shell_decay all mutations.pdf", 
                    plot_raf1, device = grDevices::cairo_pdf, height = 2.5, width = 2.5)
    plot_k13 <- multimodalallostery_plot_energy_distance_decay_expfit_contact_shell(input = "./weights_Binding_K13.txt", 
                                                                assay_sele = "K13", contact_shell = contact_shell)
    print(plot_k13)
    ggplot2::ggsave("./results/f3e_K13_contact_shell_decay all mutations.pdf", 
                    plot_k13, device = grDevices::cairo_pdf, height = 2.5, width = 2.5)
  })
  invisible(result)
}



run_figure3_e()



# ✅
# ---- Figure3F ----
run_figure3_f <- function() {
  result <- local({
    contact_shell <- data.table::fread("./5binder_contact_shell2.csv")
    plot_raf1 <- multimodalallostery_plot_energy_distance_decay_expfit_directional_no_fdr(input = "./weights_Binding_RAF.txt", 
                                                                      assay_sele = "RAF1", contact_shell = contact_shell, y_range = c(-1.5, 3))
    print(plot_raf1)
    ggplot2::ggsave("./results/f3f_RAF1_contact_shell_decay_directional.pdf", 
                    plot_raf1, device = grDevices::cairo_pdf, height = 2.5, width = 2.5)
    plot_k13 <- multimodalallostery_plot_energy_distance_decay_expfit_directional_no_fdr(input = "./weights_Binding_K13.txt", 
                                                                     assay_sele = "K13", contact_shell = contact_shell, y_range = c(-1.5, 3))
    print(plot_k13)
    ggplot2::ggsave("./results/f3f_K13_contact_shell_decay_directional.pdf", 
                    plot_k13, device = grDevices::cairo_pdf, height = 2.5, width = 2.5)
  })
  invisible(result)
}


run_figure3_f()



# ✅
# ---- Figure4B ----
run_figure4_b <- function() {
  result <- local({
    colour_scheme <- list(blue = "#1B38A6", red = "#F4270C", orange = "#F4AD0C", green = "#09B636", yellow = "#F1DD10", purple = "#C68EFD", 
                          `hot pink` = "#FF0066", `light blue` = "#75C2F6", `light red` = "#FF6A56", `dark red` = "#A31300", `dark green` = "#007A20", 
                          pink = "#FFB0A5")
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    anno <- data.table::fread("./anno_final_for_5.csv")
    K13_all_allosteric_hotspots <- c(15, 145, 10, 21, 56, 130, 139, 142, 151, 155, 178)
    K19_all_allosteric_hotspots <- c(15, 145, 10, 151, 157, 184)
    RAF1_all_allosteric_hotspots <- c(15, 16, 17, 18, 28, 32, 34, 35, 57, 60, 145, 146, 6, 10, 20, 22, 54, 55, 58, 59, 77, 144, 
                                      163, 184)
    all_allosteric_sites <- unique(c(K13_all_allosteric_hotspots, K19_all_allosteric_hotspots, RAF1_all_allosteric_hotspots))
    all_allosteric_sites <- sort(all_allosteric_sites)
    allosteric_list <- list(K13 = K13_all_allosteric_hotspots, K19 = K19_all_allosteric_hotspots, RAF1 = RAF1_all_allosteric_hotspots)
    p_triple <- multimodalallostery_plot_triple_dd_g_heatmap(ddG_K13 = "./weights_Binding_K13.txt", 
                                         ddG_K19 = "./weights_Binding_K19.txt", 
                                         ddG_RAF1 = "./weights_Binding_RAF.txt", 
                                         anno = anno, wt_aa = wt_aa, colour_scheme = colour_scheme, allosteric_sites_list = allosteric_list, legend_limits = c(-1.3, 
                                                                                                                                                               3))
    print(p_triple)
    ggplot2::ggsave("./results/f4b_triple_allosteric_sites_heatmap.pdf", 
                    p_triple, width = 10, height = 8, device = cairo_pdf)
    #ggplot2::ggsave("/triple_allosteric_sites_heatmap.png", 
    #                p_triple, width = 12, height = 10, dpi = 300)
  })
  invisible(result)
}

#install.packages("showtext")
#library(showtext)

#font_add(
#  "Arial",
#  "/System/Library/Fonts/Supplemental/Arial.ttf"
#)

#showtext_auto()
run_figure4_b()



# ✅
# ---- Figure5B ----
run_figure5_b <- function() {
  result <- local({
    COLOR_MAP <- c(`Not significant (FDR >= 0.05)` = "grey90", `Other (neutral in both)` = "grey90", `Both promoting` = "#FFB0A5", 
                   `Both disrupting` = "#F4270C", `Promoting in X / Disrupting in Y` = "#F4AD0C", `Disrupting in X / Promoting in Y` = "#F1DD10", 
                   `Allosteric only in X` = "#1B38A6", `Allosteric only in Y` = "#75C2F6")
    LEGEND_ORDER <- c("Both promoting", "Both disrupting", "Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y", 
                      "Allosteric only in X", "Allosteric only in Y", "Other (neutral in both)", "Not significant (FDR >= 0.05)")
    BINDING_SITES_MAP <- list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K55 = c(5, 24, 25, 31, 33, 36, 37, 
                                                                                                    38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74), K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 
                                                                                                                                                         71), K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138), K19 = c(68, 
                                                                                                                                                                                                                                                                                     87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137))
    NBP_RESIDUES <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    base_dir <- "./"
    anno_file <- "./anno_final_for_8.csv"
    output_dir <- "./results/"
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    cat("\n========== Analyzing RAF1 vs K13 ==========\n")
    analysis_RAF1_K13 <- multimodalallostery_analyze_protein_pair(protein_x = "RAF1", protein_y = "K13", input_file_x = file.path(base_dir, "weights_Binding_RAF.txt"), 
                                              input_file_y = file.path(base_dir, "weights_Binding_K13.txt"), anno_file = anno_file, binding_sites_map = BINDING_SITES_MAP, nbp_residues = NBP_RESIDUES, legend_order = LEGEND_ORDER, verbose = TRUE)
    plot_RAF1_K13 <- multimodalallostery_plot_protein_pair(analysis_RAF1_K13, color_map = COLOR_MAP, legend_order = LEGEND_ORDER, point_size = 2.5, alpha = 0.7, base_size = 15)
    print(plot_RAF1_K13)
    ggplot2::ggsave(file.path(output_dir, "f5b_RAF1_vs_K13_scatter_plot_NBP_triangle.pdf"), plot = plot_RAF1_K13, width = 4.5, height = 6, 
                    device = grDevices::cairo_pdf)
    cat("\n========== Analyzing RAF1 vs K55 ==========\n")
    analysis_RAF1_K55 <- multimodalallostery_analyze_protein_pair(protein_x = "RAF1", protein_y = "K55", input_file_x = file.path(base_dir, "weights_Binding_RAF.txt"), 
                                              input_file_y = file.path(base_dir, "weights_Binding_K55.txt"), anno_file = anno_file, binding_sites_map = BINDING_SITES_MAP, nbp_residues = NBP_RESIDUES, legend_order = LEGEND_ORDER, verbose = TRUE)
    plot_RAF1_K55 <- multimodalallostery_plot_protein_pair(analysis_RAF1_K55, color_map = COLOR_MAP, legend_order = LEGEND_ORDER, point_size = 2.5, alpha = 0.7, base_size = 15)
    print(plot_RAF1_K55)
    ggplot2::ggsave(file.path(output_dir, "f5b_RAF1_vs_K55_scatter_plot_NBP_triangle.pdf"), plot = plot_RAF1_K55, width = 4.5, height = 6, 
                    device = grDevices::cairo_pdf)
    cat("\n========== Analyzing RAF1 vs K27 ==========\n")
    analysis_RAF1_K27 <- multimodalallostery_analyze_protein_pair(protein_x = "RAF1", protein_y = "K27", input_file_x = file.path(base_dir, "weights_Binding_RAF.txt"), 
                                              input_file_y = file.path(base_dir, "weights_Binding_K27.txt"), anno_file = anno_file, binding_sites_map = BINDING_SITES_MAP, nbp_residues = NBP_RESIDUES, legend_order = LEGEND_ORDER, verbose = TRUE)
    plot_RAF1_K27 <- multimodalallostery_plot_protein_pair(analysis_RAF1_K27, color_map = COLOR_MAP, legend_order = LEGEND_ORDER, point_size = 2.5, alpha = 0.7, base_size = 15)
    print(plot_RAF1_K27)
    ggplot2::ggsave(file.path(output_dir, "f5b_RAF1_vs_K27_scatter_plot_NBP_triangle.pdf"), plot = plot_RAF1_K27, width = 4.5, height = 6, 
                    device = grDevices::cairo_pdf)
    cat("\n========== Analyzing K27 vs K13 ==========\n")
    analysis_K27_K13 <- multimodalallostery_analyze_protein_pair(protein_x = "K27", protein_y = "K13", input_file_x = file.path(base_dir, "weights_Binding_K27.txt"), 
                                             input_file_y = file.path(base_dir, "weights_Binding_K13.txt"), anno_file = anno_file, binding_sites_map = BINDING_SITES_MAP, nbp_residues = NBP_RESIDUES, legend_order = LEGEND_ORDER, verbose = TRUE)
    plot_K27_K13 <- multimodalallostery_plot_protein_pair(analysis_K27_K13, color_map = COLOR_MAP, legend_order = LEGEND_ORDER, point_size = 2.5, alpha = 0.7, base_size = 15)
    print(plot_K27_K13)
    ggplot2::ggsave(file.path(output_dir, "f5b_K27_vs_K13_scatter_plot_NBP_triangle.pdf"), plot = plot_K27_K13, width = 4.5, height = 6, device = grDevices::cairo_pdf)
    cat("\n========== Analyzing K13 vs K19 ==========\n")
    analysis_K13_K19 <- multimodalallostery_analyze_protein_pair(protein_x = "K13", protein_y = "K19", input_file_x = file.path(base_dir, "weights_Binding_K13.txt"), 
                                             input_file_y = file.path(base_dir, "weights_Binding_K19.txt"), anno_file = anno_file, binding_sites_map = BINDING_SITES_MAP, nbp_residues = NBP_RESIDUES, legend_order = LEGEND_ORDER, verbose = TRUE)
    plot_K13_K19 <- multimodalallostery_plot_protein_pair(analysis_K13_K19, color_map = COLOR_MAP, legend_order = LEGEND_ORDER, point_size = 2.5, alpha = 0.7, base_size = 15)
    print(plot_K13_K19)
    ggplot2::ggsave(file.path(output_dir, "f5b_K13_vs_K19_scatter_plot_NBP_triangle.pdf"), plot = plot_K13_K19, width = 4.5, height = 6, device = grDevices::cairo_pdf)
  })
  invisible(result)
}

run_figure5_b()



# ✅
# ---- Figure6B ----
run_figure6_b <- function() {
  result <- local({
    PAIR_ORDER <- c("RAF1 vs K13")
    binding_sites_map <- list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138))
    core_residues <- c(4, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 40, 42, 44, 46, 51, 52, 53, 54, 55, 
                       56, 57, 58, 68, 71, 72, 75, 77, 78, 79, 80, 81, 82, 83, 84, 89, 90, 92, 93, 96, 97, 99, 100, 101, 103, 109, 110, 111, 
                       112, 113, 114, 115, 116, 118, 125, 130, 133, 134, 137, 139, 141, 142, 143, 144, 145, 146, 151, 152, 155, 156, 157, 158, 
                       159, 160, 162, 163)
    NBP <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    functional_loop_residues <- c(10, 11, 12, 13, 14, 15, 16, 17, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
                                  40, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76)
    beta_sheets <- c(3:9, 38:44, 51:57, 77:84, 109:115, 139:143)
    surface_residues <- c(1, 2, 3, 5, 12, 13, 25:39, 41, 43, 45, 47:50, 59:67, 69, 70, 73, 74, 76, 85:88, 91, 94, 95, 98, 102, 
                          104:108, 117, 119:124, 126:129, 131, 132, 135, 136, 138, 140, 147:150, 153, 154, 161, 164:188)
    alpha_helices <- c(15:24, 67:73, 87:104, 127:136, 148:166)
    input_files <- list(RAF1 = "./weights_Binding_RAF.txt", 
                        K13 = "./weights_Binding_K13.txt")
    anno <- "./anno_final_for_5.csv"
    structure_regions <- list(Core = core_residues, Surface = surface_residues, NBP = NBP, `Functional Loop` = functional_loop_residues, 
                              `Beta Sheets` = beta_sheets, `Alpha Helices` = alpha_helices)
    target_pair <- "RAF1 vs K13"
    res <- multimodalallostery_run_pair_region_enrichment(target_pair, input_files, anno, structure_regions, binding_sites_map)
    prepared_plot <- multimodalallostery_prepare_region_enrichment_plot_data(res, target_pair, names(structure_regions))
    color_map <- c(Correlated = "#F4AD0C", `Anti-correlated` = "#1B38A6", Other = "grey80")
    p <- multimodalallostery_plot_region_enrichment(prepared_plot, color_map)
    print(p)
    output_dir <- "./results/"
    ggplot2::ggsave(paste0(output_dir, "f6b_barplot_RAF1_vs_K13_core_surface_nbp_fl_bsheet_ahelix_enrichment.pdf"), plot = p, width = 10, 
                    height = 6, device = grDevices::cairo_pdf)
    outputs <- multimodalallostery_format_region_enrichment_outputs(res)
    data.table::fwrite(outputs$odds, paste0(output_dir, "f6b_OR_values_RAF1_vs_K13_all_six_regions.csv"))
    data.table::fwrite(outputs$fractions, paste0(output_dir, "f6b_fraction_summary_RAF1_vs_K13_all_six_regions.csv"))
    cat("\n\nAnalysis completed for", target_pair, "\n")
    cat("Results saved to:", output_dir, "\n")
    cat("\nAnalyzed regions:\n")
    for (region in names(structure_regions)) {
      cat("  -", region, ":", length(structure_regions[[region]]), "residues\n")
    }
  })
  invisible(result)
}


run_figure6_b()


# ✅
# ---- FigureS1B ----
run_figure_s1_b <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    result_abundance <- multimodalallostery_compare_fitness_libraries_singlemut_overall_no_block1(lib1_block2 = "./fitness_RData/CW_RAS_abundance_2_fitness_replicates_fullseq.RData", 
                                                                              lib1_block3 = "./fitness_RData/CW_RAS_abundance_3_fitness_replicates_fullseq.RData", 
                                                                              lib2_block2 = "./fitness_RData/MA_RAS_abundance_2_fitness_replicates_fullseq.RData", 
                                                                              lib2_block3 = "./fitness_RData/MA_RAS_abundance_3_fitness_replicates_fullseq.RData", 
                                                                              wt_aa = wt_aa, output_file = "./results/s1b_comparison_of_fitness_data_Abundance_Overall_no_block1.pdf", 
                                                                              x_lab = "Abundance nicking library fitness", y_lab = "Abundance synthetic library fitness", main_title = "Comparison of fitness data between synthetic library and nicking library\nsingle mutation", 
                                                                              point_alpha = 0.3, plot_width = 5, plot_height = 5)
    result_raf1 <- multimodalallostery_compare_fitness_libraries_singlemut_overall_no_block1(lib1_block2 = "./fitness_RData/CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData", 
                                                                         lib1_block3 = "./fitness_RData/CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData", 
                                                                         lib2_block2 = "./fitness_RData/MA_RAS_binding_RAF_2_fitness_replicates_fullseq.RData", 
                                                                         lib2_block3 = "./fitness_RData/MA_RAS_binding_RAF_3_fitness_replicates_fullseq.RData", 
                                                                         wt_aa = wt_aa, output_file = "./results/s1b_comparison_of_fitness_data_RAF1_Overall_no_block1.pdf", 
                                                                         x_lab = "RAF1 nicking library fitness", y_lab = "RAF1 synthetic library fitness", main_title = "Comparison of fitness data between synthetic library and nicking library\nsingle mutation", 
                                                                         point_alpha = 0.3, plot_width = 5, plot_height = 5)
    result_K55 <- multimodalallostery_compare_fitness_libraries_singlemut_overall_no_block1(lib1_block2 = "./fitness_RData/CW_RAS_binding_K55_2_fitness_replicates_fullseq.RData", 
                                                                        lib1_block3 = "./fitness_RData/CW_RAS_binding_K55_3_fitness_replicates_fullseq.RData", 
                                                                        lib2_block2 = "./fitness_RData/MA_RAS_binding_K55_2_fitness_replicates_fullseq.RData", 
                                                                        lib2_block3 = "./fitness_RData/MA_RAS_binding_K55_3_fitness_replicates_fullseq.RData", 
                                                                        wt_aa = wt_aa, output_file = "./results/s1b_comparison_of_fitness_data_K55_Overall_no_block1.pdf", 
                                                                        x_lab = "K55 nicking library fitness", y_lab = "K55 synthetic library fitness", main_title = "Comparison of fitness data between synthetic library and nicking library\nsingle mutation", 
                                                                        point_alpha = 0.3, plot_width = 5, plot_height = 5)
    result_K27 <- multimodalallostery_compare_fitness_libraries_singlemut_overall_no_block1(lib1_block2 = "./fitness_RData/CW_RAS_binding_K27_2_fitness_replicates_fullseq.RData", 
                                                                        lib1_block3 = "./fitness_RData/CW_RAS_binding_K27_3_fitness_replicates_fullseq.RData", 
                                                                        lib2_block2 = "./fitness_RData/MA_RAS_binding_K27_2_fitness_replicates_fullseq.RData", 
                                                                        lib2_block3 = "./fitness_RData/MA_RAS_binding_K27_3_fitness_replicates_fullseq.RData", 
                                                                        wt_aa = wt_aa, output_file = "./results/s1b_comparison_of_fitness_data_K27_Overall_no_block1.pdf", 
                                                                        x_lab = "K27 nicking library fitness", y_lab = "K27 synthetic library fitness", main_title = "Comparison of fitness data between synthetic library and nicking library\nsingle mutation", 
                                                                        point_alpha = 0.3, plot_width = 5, plot_height = 5)
  })
  invisible(result)
}


run_figure_s1_b()


# ✅
# ---- FigureS1C ----
run_figure_s1_c <- function() {
  result <- local({
    colour_scheme <- list(blue = "#1B38A6", red = "#F4270C", orange = "#F4AD0C", green = "#09B636", yellow = "#F1DD10", purple = "#C68EFD", 
                          `hot pink` = "#FF0066", `light blue` = "#75C2F6", `light red` = "#FF6A56", `dark red` = "#A31300", `dark green` = "#007A20", 
                          pink = "#FFB0A5")
    K13_nor_df <- krasddpcams::krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./fitness_RData/MA_RAS_binding_K13_1_fitness_replicates_fullseq.RData", 
                                                                         block2_dimsum_df = "./fitness_RData/MA_RAS_binding_K13_2_fitness_replicates_fullseq.RData", 
                                                                         block3_dimsum_df = "./fitness_RData/MA_RAS_binding_K13_3_fitness_replicates_fullseq.RData")
    K19_nor_df <- krasddpcams::krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./fitness_RData/MA_RAS_binding_K19_1_fitness_replicates_fullseq.RData", 
                                                                         block2_dimsum_df = "./fitness_RData/MA_RAS_binding_K19_2_fitness_replicates_fullseq.RData", 
                                                                         block3_dimsum_df = "./fitness_RData/MA_RAS_binding_K19_3_fitness_replicates_fullseq.RData")
    stability_nor_df <- krasddpcams::krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./fitness_RData/CW_RAS_abundance_1_fitness_replicates_fullseq.RData", 
                                                                               block2_dimsum_df = "./fitness_RData/MA_RAS_abundance_2_fitness_replicates_fullseq.RData", 
                                                                               block3_dimsum_df = "./fitness_RData/MA_RAS_abundance_3_fitness_replicates_fullseq.RData")
    RAF1_nor_df <- krasddpcams::krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./fitness_RData/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData", 
                                                                          block2_dimsum_df = "./fitness_RData/MA_RAS_binding_RAF_2_fitness_replicates_fullseq.RData", 
                                                                          block3_dimsum_df = "./fitness_RData/MA_RAS_binding_RAF_3_fitness_replicates_fullseq.RData")
    K55_nor_df <- krasddpcams::krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./fitness_RData/CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData", 
                                                                         block2_dimsum_df = "./fitness_RData/MA_RAS_binding_K55_2_fitness_replicates_fullseq.RData", 
                                                                         block3_dimsum_df = "./fitness_RData/MA_RAS_binding_K55_3_fitness_replicates_fullseq.RData")
    K27_nor_df <- krasddpcams::krasddpcams__normalize_growthrate_fitness(block1_dimsum_df = "./fitness_RData/CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData", 
                                                                         block2_dimsum_df = "./fitness_RData/MA_RAS_binding_K27_2_fitness_replicates_fullseq.RData", 
                                                                         block3_dimsum_df = "./fitness_RData/MA_RAS_binding_K27_3_fitness_replicates_fullseq.RData")
    d_K13 <- multimodalallostery_plot_fitness_correlation_blocks_BI2(K13_nor_df, "K13", colour_scheme)
    print(d_K13)
    ggplot2::ggsave("./results/s1c_fitness_correlation_blocks_K13.pdf", 
                    d_K13, device = grDevices::cairo_pdf, height = 4, width = 4)
    d_K19 <- multimodalallostery_plot_fitness_correlation_blocks_BI2(K19_nor_df, "K19", colour_scheme)
    print(d_K19)
    ggplot2::ggsave("./results/s1c_fitness_correlation_blocks_K19.pdf", 
                    d_K19, device = grDevices::cairo_pdf, height = 4, width = 4)
    d_stability <- multimodalallostery_plot_fitness_correlation_blocks_BI1(stability_nor_df, "stability", colour_scheme)
    print(d_stability)
    ggplot2::ggsave("./results/s1c_fitness_correlation_blocks_stability.pdf", 
                    d_stability, device = grDevices::cairo_pdf, height = 4, width = 4)
    d_RAF1 <- multimodalallostery_plot_fitness_correlation_blocks_BI1(RAF1_nor_df, "RAF1", colour_scheme)
    print(d_RAF1)
    ggplot2::ggsave("./results/s1c_fitness_correlation_blocks_RAF1.pdf", 
                    d_RAF1, device = grDevices::cairo_pdf, height = 4, width = 4)
    d_K55 <- multimodalallostery_plot_fitness_correlation_blocks_BI1(K55_nor_df, "K55", colour_scheme)
    print(d_K55)
    ggplot2::ggsave("./results/s1c_fitness_correlation_blocks_K55.pdf", 
                    d_K55, device = grDevices::cairo_pdf, height = 4, width = 4)
    d_K27 <- multimodalallostery_plot_fitness_correlation_blocks_BI1(K27_nor_df, "K27", colour_scheme)
    print(d_K27)
    ggplot2::ggsave("./results/s1c_fitness_correlation_blocks_K27.pdf", 
                    d_K27, device = grDevices::cairo_pdf, height = 4, width = 4)
  })
  invisible(result)
}


run_figure_s1_c()



# ✅
# ---- FigureS1D ----
run_figure_s1_d <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    nor_fit <- wlab.block::nor_fitness(block1 = "./fitness_RData/CW_RAS_abundance_1_fitness_replicates_fullseq.RData", 
                                       block2 = "./fitness_RData/MA_RAS_abundance_2_fitness_replicates_fullseq.RData", 
                                       block3 = "./fitness_RData/MA_RAS_abundance_3_fitness_replicates_fullseq.RData")
    nor_fit_single <- wlab.block::nor_fitness_single_mut(input = nor_fit)
    nor_fit_single <- wlab.block::pos_id(nor_fit_single, wt_aa)
    multimodalallostery_fitness_heatmap(nor_fit_single, wt_aa, title = "KRAS-Abundance", legend_limits = c(-2, 1.5))
    ggplot2::ggsave("./results/s1d_KRAS-Abundance.pdf", 
                    height = 6, width = 20)
    nor_fit <- wlab.block::nor_fitness(block1 = "./fitness_RData/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData", 
                                       block2 = "./fitness_RData/MA_RAS_binding_RAF_2_fitness_replicates_fullseq.RData", 
                                       block3 = "./fitness_RData/MA_RAS_binding_RAF_3_fitness_replicates_fullseq.RData")
    nor_fit_single <- wlab.block::nor_fitness_single_mut(input = nor_fit)
    nor_fit_single <- wlab.block::pos_id(nor_fit_single, wt_aa)
    multimodalallostery_fitness_heatmap(nor_fit_single, wt_aa, title = "KRAS-RAF1", legend_limits = c(-2, 1.5))
    ggplot2::ggsave("./results/s1d_KRAS-RAF1.pdf", height = 6, 
                    width = 20)
    nor_fit <- wlab.block::nor_fitness(block1 = "./fitness_RData/MA_RAS_binding_K13_1_fitness_replicates_fullseq.RData", 
                                       block2 = "./fitness_RData/MA_RAS_binding_K13_2_fitness_replicates_fullseq.RData", 
                                       block3 = "./fitness_RData/MA_RAS_binding_K13_3_fitness_replicates_fullseq.RData")
    nor_fit_single <- wlab.block::nor_fitness_single_mut(input = nor_fit)
    nor_fit_single <- wlab.block::pos_id(nor_fit_single, wt_aa)
    multimodalallostery_fitness_heatmap(nor_fit_single, wt_aa, title = "KRAS-DARPin K13", legend_limits = c(-2, 1.5))
    ggplot2::ggsave("./results/s1d_KRAS-DARPin K13.pdf", 
                    height = 6, width = 20)
    nor_fit <- wlab.block::nor_fitness(block1 = "./fitness_RData/MA_RAS_binding_K19_1_fitness_replicates_fullseq.RData", 
                                       block2 = "./fitness_RData/MA_RAS_binding_K19_2_fitness_replicates_fullseq.RData", 
                                       block3 = "./fitness_RData/MA_RAS_binding_K19_3_fitness_replicates_fullseq.RData")
    nor_fit_single <- wlab.block::nor_fitness_single_mut(input = nor_fit)
    nor_fit_single <- wlab.block::pos_id(nor_fit_single, wt_aa)
    multimodalallostery_fitness_heatmap(nor_fit_single, wt_aa, title = "KRAS-DARPin K19", legend_limits = c(-2, 1.5))
    ggplot2::ggsave("./results/s1d_KRAS-DARPin K19.pdf", 
                    height = 6, width = 20)
    nor_fit <- wlab.block::nor_fitness(block1 = "./fitness_RData/CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData", 
                                       block2 = "./fitness_RData/MA_RAS_binding_K55_2_fitness_replicates_fullseq.RData", 
                                       block3 = "./fitness_RData/MA_RAS_binding_K55_3_fitness_replicates_fullseq.RData")
    nor_fit_single <- wlab.block::nor_fitness_single_mut(input = nor_fit)
    nor_fit_single <- wlab.block::pos_id(nor_fit_single, wt_aa)
    multimodalallostery_fitness_heatmap(nor_fit_single, wt_aa, title = "KRAS-DARPin K55", legend_limits = c(-2, 1.5))
    ggplot2::ggsave("./results/s1d_KRAS-DARPin K55.pdf", 
                    height = 6, width = 20)
    nor_fit <- wlab.block::nor_fitness(block1 = "./fitness_RData/CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData", 
                                       block2 = "./fitness_RData/MA_RAS_binding_K27_2_fitness_replicates_fullseq.RData", 
                                       block3 = "./fitness_RData/MA_RAS_binding_K27_3_fitness_replicates_fullseq.RData")
    nor_fit_single <- wlab.block::nor_fitness_single_mut(input = nor_fit)
    nor_fit_single <- wlab.block::pos_id(nor_fit_single, wt_aa)
    multimodalallostery_fitness_heatmap(nor_fit_single, wt_aa, title = "KRAS-DARPin K27", legend_limits = c(-2, 1.5))
    ggplot2::ggsave("./results/s1d_KRAS-DARPin K27.pdf", 
                    height = 6, width = 20)
  })
  invisible(result)
}


run_figure_s1_d()



# ✅
# ---- FigureS2A ----
run_figure_s2_a <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    Folding_pre_ob_fitness <- multimodalallostery_merge_ddgf_fitness_blocks(prediction = "./predicted_phenotypes_all.txt", 
                                                        folding_ddG = "./weights_Folding.txt", 
                                                        block1_dimsum_df = "./fitness_RData/CW_RAS_abundance_1_fitness_replicates_fullseq.RData", 
                                                        block2_dimsum_df = "./fitness_RData/CW_RAS_abundance_2_fitness_replicates_fullseq.RData", 
                                                        block3_dimsum_df = "./fitness_RData/MA_RAS_abundance_2_fitness_replicates_fullseq.RData", 
                                                        block4_dimsum_df = "./fitness_RData/CW_RAS_abundance_3_fitness_replicates_fullseq.RData", 
                                                        block5_dimsum_df = "./fitness_RData/MA_RAS_abundance_3_fitness_replicates_fullseq.RData", 
                                                        wt_aa_input = wt_aa)
    krasddpcams::krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness, fold_n = 1, mochi_parameters = "./linears_weights_Abundance1.txt", 
                                                  phenotypen = 1, RT = 0.001987 * (273 + 30), bin_input = 50)
    ggplot2::ggsave("./results/s2a_block1_ddGf_fitness.pdf", 
                    device = grDevices::cairo_pdf, height = 35, width = 60, units = "mm")
    krasddpcams::krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness, fold_n = 1, mochi_parameters = "./linears_weights_Abundance2_1.txt", 
                                                  phenotypen = 2, RT = 0.001987 * (273 + 30), bin_input = 50)
    ggplot2::ggsave("./results/s2a_block2_1_ddGf_fitness.pdf", 
                    device = grDevices::cairo_pdf, height = 35, width = 60, units = "mm")
    krasddpcams::krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness, fold_n = 1, mochi_parameters = "./linears_weights_Abundance2_2.txt", 
                                                  phenotypen = 3, RT = 0.001987 * (273 + 30), bin_input = 50)
    ggplot2::ggsave("./results/s2a_block2_2_ddGf_fitness.pdf", 
                    device = grDevices::cairo_pdf, height = 35, width = 60, units = "mm")
    krasddpcams::krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness, fold_n = 1, mochi_parameters = "./linears_weights_Abundance3_1.txt", 
                                                  phenotypen = 4, RT = 0.001987 * (273 + 30), bin_input = 50)
    ggplot2::ggsave("./results/s2a_block3_1_ddGf_fitness.pdf", 
                    device = grDevices::cairo_pdf, height = 35, width = 60, units = "mm")
    krasddpcams::krasddpcams__plot2d_ddGf_fitness(pre_nor = Folding_pre_ob_fitness, fold_n = 1, mochi_parameters = "./linears_weights_Abundance3_2.txt", 
                                                  phenotypen = 5, RT = 0.001987 * (273 + 30), bin_input = 50)
    ggplot2::ggsave("./results/s2a_block3_2_ddGf_fitness.pdf", 
                    device = grDevices::cairo_pdf, height = 35, width = 60, units = "mm")
  })
  invisible(result)
}


run_figure_s2_a()



# ✅
# ---- FigureS2B ----
run_figure_s2_b <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    abundance_fitness_data <- multimodalallostery_merge_ddgf_fitness_blocks(prediction = "./predicted_phenotypes_all.txt", 
                                                                            folding_ddG = "./weights_Folding.txt", 
                                                                            block1_dimsum_df = "./fitness_RData/CW_RAS_abundance_1_fitness_replicates_fullseq.RData", 
                                                                            block2_dimsum_df = "./fitness_RData/CW_RAS_abundance_2_fitness_replicates_fullseq.RData", 
                                                                            block3_dimsum_df = "./fitness_RData/MA_RAS_abundance_2_fitness_replicates_fullseq.RData", 
                                                                            block4_dimsum_df = "./fitness_RData/CW_RAS_abundance_3_fitness_replicates_fullseq.RData", 
                                                                            block5_dimsum_df = "./fitness_RData/MA_RAS_abundance_3_fitness_replicates_fullseq.RData", 
                                                                            wt_aa_input = wt_aa)
    multimodalallostery_plot_ddg_fitness_per_block(pre_nor = abundance_fitness_data, phenotypen = 1, rotate_x_axis = TRUE)
    ggplot2::ggsave("./results/s2b_abundance_block1_evaluation_fitness_pre_vs_ob.pdf", 
                    device = grDevices::cairo_pdf, height = 45, width = 60, units = "mm")
    multimodalallostery_plot_ddg_fitness_per_block(pre_nor = abundance_fitness_data, phenotypen = 2, rotate_x_axis = TRUE)
    ggplot2::ggsave("./results/s2b_abundance_block2_1_evaluation_fitness_pre_vs_ob.pdf", 
                    device = grDevices::cairo_pdf, height = 45, width = 60, units = "mm")
    multimodalallostery_plot_ddg_fitness_per_block(pre_nor = abundance_fitness_data, phenotypen = 3, rotate_x_axis = TRUE)
    ggplot2::ggsave("./results/s2b_abundance_block2_2_evaluation_fitness_pre_vs_ob.pdf", 
                    device = grDevices::cairo_pdf, height = 45, width = 60, units = "mm")
    multimodalallostery_plot_ddg_fitness_per_block(pre_nor = abundance_fitness_data, phenotypen = 4, rotate_x_axis = TRUE)
    ggplot2::ggsave("./results/s2b_abundance_block3_1_evaluation_fitness_pre_vs_ob.pdf", 
                    device = grDevices::cairo_pdf, height = 45, width = 60, units = "mm")
    multimodalallostery_plot_ddg_fitness_per_block(pre_nor = abundance_fitness_data, phenotypen = 5, rotate_x_axis = TRUE)
    ggplot2::ggsave("./results/s2b_abundance_block3_2_evaluation_fitness_pre_vs_ob.pdf", 
                    device = grDevices::cairo_pdf, height = 45, width = 60, units = "mm")
  })
  invisible(result)
}


run_figure_s2_b()


# ✅
# ---- FigureS2C ----
run_figure_s2_c <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    K13_pre_ob_fitness <- multimodalallostery_k13_19_get_ob_pre_fitness_binding_correlation_3blocks(prediction = "./predicted_phenotypes_all.txt", 
                                                                                block1_dimsum_df = "./fitness_RData/MA_RAS_binding_K13_1_fitness_replicates_fullseq.RData", 
                                                                                block2_dimsum_df = "./fitness_RData/MA_RAS_binding_K13_2_fitness_replicates_fullseq.RData", 
                                                                                block3_dimsum_df = "./fitness_RData/MA_RAS_binding_K13_3_fitness_replicates_fullseq.RData", 
                                                                                assay_sele = "K13", wt_aa_input = wt_aa)
    multimodalallostery_plot_ddg_fitness_per_block(pre_nor = K13_pre_ob_fitness, phenotypen = 6, rotate_x_axis = TRUE)
    ggplot2::ggsave("./results/s2c_K13_fitness_pre_vs_ob_block1.pdf", 
                    device = grDevices::cairo_pdf, height = 45, width = 60, units = "mm")
    multimodalallostery_plot_ddg_fitness_per_block(pre_nor = K13_pre_ob_fitness, phenotypen = 7, rotate_x_axis = TRUE)
    ggplot2::ggsave("./results/s2c_K13_fitness_pre_vs_ob_block2.pdf", 
                    device = grDevices::cairo_pdf, height = 45, width = 60, units = "mm")
    multimodalallostery_plot_ddg_fitness_per_block(pre_nor = K13_pre_ob_fitness, phenotypen = 8, rotate_x_axis = TRUE)
    ggplot2::ggsave("./results/s2c_K13_fitness_pre_vs_ob_block3.pdf", 
                    device = grDevices::cairo_pdf, height = 45, width = 60, units = "mm")
  })
  invisible(result)
}


run_figure_s2_c()

# ✅
# ---- FigureS2D ----
run_figure_s2_d <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    K19_pre_ob_fitness <- multimodalallostery_k13_19_get_ob_pre_fitness_binding_correlation_3blocks(prediction = "./predicted_phenotypes_all.txt", 
                                                                                                    block1_dimsum_df = "./fitness_RData/MA_RAS_binding_K19_1_fitness_replicates_fullseq.RData", 
                                                                                                    block2_dimsum_df = "./fitness_RData/MA_RAS_binding_K19_2_fitness_replicates_fullseq.RData", 
                                                                                                    block3_dimsum_df = "./fitness_RData/MA_RAS_binding_K19_3_fitness_replicates_fullseq.RData", 
                                                                                assay_sele = "K19", wt_aa_input = wt_aa)
    multimodalallostery_plot_ddg_fitness_per_block(pre_nor = K19_pre_ob_fitness, phenotypen = 9, rotate_x_axis = TRUE)
    ggplot2::ggsave("./results/s2d_K19_fitness_pre_vs_ob_block1.pdf", 
                    device = grDevices::cairo_pdf, height = 45, width = 60, units = "mm")
    multimodalallostery_plot_ddg_fitness_per_block(pre_nor = K19_pre_ob_fitness, phenotypen = 10, rotate_x_axis = TRUE)
    ggplot2::ggsave("./results/s2d_K19_fitness_pre_vs_ob_block2.pdf", 
                    device = grDevices::cairo_pdf, height = 45, width = 60, units = "mm")
    multimodalallostery_plot_ddg_fitness_per_block(pre_nor = K19_pre_ob_fitness, phenotypen = 11, rotate_x_axis = TRUE)
    ggplot2::ggsave("./results/s2d_K19_fitness_pre_vs_ob_block3.pdf", 
                    device = grDevices::cairo_pdf, height = 45, width = 60, units = "mm")
  })
  invisible(result)
}


run_figure_s2_d()


# ✅
# ---- FigureS2F ----
run_figure_s2_f <- function() {
  result <- local({
    ddG_folding_new <- data.table::fread("./weights_Folding.txt")
    ddG_folding_new <- ddG_folding_new[, c(1, 3, 20:22)]
    ddG_folding_weng <- data.table::fread("./weights_Folding_weng.txt")
    ddG_folding_weng <- ddG_folding_weng[, c(1, 3, 20:22)]
    data1 = ddG_folding_new
    data2 = ddG_folding_weng
    correlation_plot <- multimodalallostery_plot_dd_g_correlation(data1 = ddG_folding_new, data2 = ddG_folding_weng, data1_name = "Folding_this study", 
                                              data2_name = "Folding_Weng", output_file = "./results/s2f_ddG correlation between folding_new and folding_weng.pdf", 
                                              limits = c(-1.6, 2.8))
    print(correlation_plot)
    ddG_RAF1_new <- data.table::fread("./weights_Binding_RAF.txt")
    ddG_RAF1_new <- ddG_RAF1_new[, c(1, 3, 20:22)]
    ddG_RAF1_weng <- data.table::fread("./weights_Binding_RAF_weng.txt")
    ddG_RAF1_weng <- ddG_RAF1_weng[, c(1, 3, 20:22)]
    data1 = ddG_RAF1_new
    data2 = ddG_RAF1_weng
    correlation_plot <- multimodalallostery_plot_dd_g_correlation(data1 = ddG_RAF1_new, data2 = ddG_RAF1_weng, data1_name = "RAF1_this study", data2_name = "RAF1_Weng", 
                                              output_file = "./results/s2f_ddG correlation between RAF1_new and folding_weng.pdf", 
                                              limits = c(-1.6, 2.8))
    print(correlation_plot)
    ddG_K55_new <- data.table::fread("./weights_Binding_K55.txt")
    ddG_K55_new <- ddG_K55_new[, c(1, 3, 20:22)]
    ddG_K55_weng <- data.table::fread("./weights_Binding_K55_weng.txt")
    ddG_K55_weng <- ddG_K55_weng[, c(1, 3, 20:22)]
    data1 = ddG_K55_new
    data2 = ddG_K55_weng
    correlation_plot <- multimodalallostery_plot_dd_g_correlation(data1 = ddG_K55_new, data2 = ddG_K55_weng, data1_name = "K55_this study", data2_name = "K55_Weng", 
                                              output_file = "./results/s2f_ddG correlation between K55_new and folding_weng.pdf", 
                                              limits = c(-1.6, 2.8))
    print(correlation_plot)
    ddG_K27_new <- data.table::fread("./weights_Binding_K27.txt")
    ddG_K27_new <- ddG_K27_new[, c(1, 3, 20:22)]
    ddG_K27_weng <- data.table::fread("./weights_Binding_K27_weng.txt")
    ddG_K27_weng <- ddG_K27_weng[, c(1, 3, 20:22)]
    data1 = ddG_K27_new
    data2 = ddG_K27_weng
    correlation_plot <- multimodalallostery_plot_dd_g_correlation(data1 = ddG_K27_new, data2 = ddG_K27_weng, data1_name = "K27_this study", data2_name = "K27_Weng", 
                                              output_file = "./results/s2f_ddG correlation between K27_new and folding_weng.pdf", 
                                              limits = c(-1.6, 2.8))
    print(correlation_plot)
  })
  invisible(result)
}


run_figure_s2_f()




# ✅
# ---- FigureS3A ----
run_figure_s3_a <- function() {
  result <- local({
    annotation_file <- "./anno_final_for_8.csv"
    new_ddg_file <- "./weights_Binding_RAF.txt"
    beta_sheet_ranges <- data.frame(xstart = c(3, 38, 51, 77, 109, 139), xend = c(9, 44, 57, 84, 115, 143), col = paste0("b", 
                                                                                                                         1:6))
    annotation <- data.table::fread(annotation_file)
    new_plot <- multimodalallostery_plot_beta_sheet_ddg(multimodalallostery_prepare_beta_sheet_ddg(data.table::fread(new_ddg_file), annotation, beta_sheet_ranges), "Binding free energy change(RAF1) \n(kcal/mol)")
    new_plot
    ggplot2::ggsave("./results/s3a_ddG_betasheet decay new RAF1 energy data.pdf", 
                    plot = new_plot, device = grDevices::cairo_pdf, height = 4, width = 4)
  })
  invisible(result)
}


run_figure_s3_a()




# ✅
# ---- FigureS3B ----
run_figure_s3_b <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    rects_sheet <- data.frame(xstart = c(3, 38, 51, 77, 109, 139), xend = c(9, 44, 57, 84, 115, 143), col = c("β1", "β2", "β3", 
                                                                                                              "β4", "β5", "β6"))
    rects_helix <- data.frame(xstart = c(15, 67, 87, 127, 148), xend = c(24, 73, 104, 136, 166), col = c("α1", "α2", "α3", 
                                                                                                         "α4", "α5"))
    multimodalallostery_manhatta_plot_single_assay(input = "./weights_Binding_K19.txt", 
                               assay_sele = "K19", anno = "./anno_final_for_5.csv", 
                               rects_sheet = rects_sheet, rects_alpha = rects_helix, wt_aa = wt_aa)
    ggplot2::ggsave("./results/s3b_K19_Manhattan_plot.pdf", 
                    device = grDevices::cairo_pdf, height = 9, width = 12)
  })
  invisible(result)
}



run_figure_s3_b()



# ✅
# ---- FigureS3C ----
run_figure_s3_c <- function() {
  result <- local({
    wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
    NBP_positions <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    beta_sheets <- c(3:9, 38:44, 51:57, 77:84, 109:115, 139:143)
    contact_shell_file <- "./5binder_contact_shell2.csv"
    output_dir <- "./results/"
    result_K13 <- multimodalallostery_process_single_assay(input = "./weights_Binding_K13.txt", 
                                       assay_sele = "K13", anno = "./anno_final_for_5.csv", 
                                       contact_shell_file = contact_shell_file, wt_aa = wt_aa, NBP_positions = NBP_positions, beta_sheets = beta_sheets)
    result_K19 <- multimodalallostery_process_single_assay(input = "./weights_Binding_K19.txt", 
                                       assay_sele = "K19", anno = "./anno_final_for_5.csv", 
                                       contact_shell_file = contact_shell_file, wt_aa = wt_aa, NBP_positions = NBP_positions, beta_sheets = beta_sheets)
    result_RAF1 <- multimodalallostery_process_single_assay(input = "./weights_Binding_RAF.txt", 
                                        assay_sele = "RAF1", anno = "./anno_final_for_5.csv", 
                                        contact_shell_file = contact_shell_file, wt_aa = wt_aa, NBP_positions = NBP_positions, beta_sheets = beta_sheets)
    plot_data_K13 <- multimodalallostery_prepare_plot_data(result_K13$NBP_results, result_K13$second_shell_results, result_K13$beta_sheet_results, 
                                       "K13")
    plot_K13 <- multimodalallostery_create_enrichment_plot(plot_data_K13, "K13", output_dir)
    plot_K13
    plot_data_K19 <- multimodalallostery_prepare_plot_data(result_K19$NBP_results, result_K19$second_shell_results, result_K19$beta_sheet_results, 
                                       "K19")
    plot_K19 <- multimodalallostery_create_enrichment_plot(plot_data_K19, "K19", output_dir)
    plot_K19
    plot_data_RAF1 <- multimodalallostery_prepare_plot_data(result_RAF1$NBP_results, result_RAF1$second_shell_results, result_RAF1$beta_sheet_results, 
                                        "RAF1")
    plot_RAF1 <- multimodalallostery_create_enrichment_plot(plot_data_RAF1, "RAF1", output_dir)
    plot_RAF1
    cat("\n\n", paste0(rep("=", 80), collapse = ""), "\n")
    cat("PLOT DATA SUMMARY\n")
    cat(paste0(rep("=", 80), collapse = ""), "\n\n")
    cat("K13 Plot Data:\n")
    print(plot_data_K13)
    cat("\nK19 Plot Data:\n")
    print(plot_data_K19)
    cat("\nRAF1 Plot Data:\n")
    print(plot_data_RAF1)
    cat("\n\nPlots saved to:\n")
    cat("  -", file.path(output_dir, "s3c_Enrichment_plot_K13.pdf\n"))
    cat("  -", file.path(output_dir, "s3c_Enrichment_plot_K19.pdf\n"))
    cat("  -", file.path(output_dir, "s3c_Enrichment_plot_RAF1.pdf\n"))
  })
  invisible(result)
}



run_figure_s3_c()


# ✅
# ---- FigureS3E ----
run_figure_s3_e <- function() {
  result <- local({
    plot_k19 <- multimodalallostery_plot_energy_distance_decay_expfit_annotation(input = "./weights_Binding_K19.txt", 
                                                             assay_sele = "K19", anno_file = "./anno_final_for_5.csv")
    print(plot_k19)
    ggplot2::ggsave("./results/s3e_K19_distance_decay all mutations.pdf", 
                    plot_k19, device = grDevices::cairo_pdf, width = 2.5, height = 2.5, units = "in", dpi = 300)
  })
  invisible(result)
}

run_figure_s3_e()




# ---- FigureS3F ----
run_figure_s3_f <- function() {
  result <- local({
    plot_K19 <- multimodalallostery_plot_energy_distance_decay_directional_no_filter(input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt", 
                                                                 assay_sele = "K19", anno_file = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv")
    print(plot_K19)
    ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/test_20260810/K19_directional_decay.pdf", 
                    plot_K19, device = grDevices::cairo_pdf, width = 2.5, height = 2.5, units = "in", dpi = 300)
  })
  invisible(result)
}


run_figure_s3_f()



# ✅
# ---- FigureS3G ----
run_figure_s3_g <- function() {
  result <- local({
    contact_shell <- data.table::fread("./5binder_contact_shell2.csv")
    plot_k19 <- multimodalallostery_plot_energy_distance_decay_expfit_contact_shell(input = "./weights_Binding_K19.txt", 
                                                                assay_sele = "K19", contact_shell = contact_shell)
    print(plot_k19)
    ggplot2::ggsave("./results/s3g_K19_contact_shell_decay all mutations.pdf", 
                    plot_k19, device = grDevices::cairo_pdf, height = 2.5, width = 2.5)
  })
  invisible(result)
}


run_figure_s3_g()



# ✅
# ---- FigureS3H ----
run_figure_s3_h <- function() {
  result <- local({
    contact_shell <- data.table::fread("./5binder_contact_shell2.csv")
    plot_k19 <- multimodalallostery_plot_energy_distance_decay_expfit_directional_no_fdr(input = "./weights_Binding_K19.txt", 
                                                                     assay_sele = "K19", contact_shell = contact_shell, y_range = c(-1.5, 3))
    print(plot_k19)
    ggplot2::ggsave("./results/s3h_K19_contact_shell_decay_directional.pdf", 
                    plot_k19, device = grDevices::cairo_pdf, height = 2.5, width = 2.5)
  })
  invisible(result)
}

run_figure_s3_h()


# ✅
# ---- FigureS4A ----
run_figure_s4_a <- function() {
  result <- local({
    ddG_files <- c("./weights_Binding_K13.txt", 
                   "./weights_Binding_K19.txt", 
                   "./weights_Binding_RAF.txt")
    assays <- c("K13", "K19", "RAF1")
    anno_file <- "./anno_final_for_5.csv"
    result <- multimodalallostery_plot_weighted_mean_dd_g_distance_with_cross(ddG_files = ddG_files, assays = assays, anno_file = anno_file, x_intercept = 5, 
                                                          output_file = "./results/s4a_K13_K19_RAF1_with_cross_labels_colored_text.pdf", 
                                                          base_font_size = 15, point_size = 0.8, text_repel_size = 5)
    #print(result$plot)
    print(paste("Threshold value:", result$threshold))
    cat("\n\nAdditionally annotated cross-hotspots (derived from other proteins, rather than intrinsic allosteric hotspots):\n")
    print(result$cross_hotspots[, .(assay, Pos_real, cross_hotspot_source, mean, distance_bp, site_type)])
    data.table::fwrite(result$cross_hotspots[, .(assay, Pos_real, cross_hotspot_source, mean, sigma, distance_bp, site_type, count)], "./results/s4a_cross_hotspots_data.csv")
  })
  invisible(result)
}


run_figure_s4_a()



# ✅
# ---- FigureS5A ----
run_figure_s5_a <- function() {
  result <- local({
    binding_sites_map <- list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), RALGDS = c(24, 25, 31, 33, 36, 37, 
                                                                                                       38, 39, 40, 41, 56, 64, 67), PI3KCG = c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73), SOS1 = c(1, 22, 24, 
                                                                                                                                                                                                                    25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153), K55 = c(5, 24, 25, 
                                                                                                                                                                                                                                                                                                                               31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74), K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 
                                                                                                                                                                                                                                                                                                                                                                                                    52, 54, 67, 70, 71), K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 
                                                                                                                                                                                                                                                                                                                                                                                                                                 137, 138), K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137))
    BINDER_ORDER <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")
    input_files <- list(RAF1 = "./weights_Binding_RAF.txt", 
                        RALGDS = "./weights_Binding_RAL.txt", 
                        PI3KCG = "./weights_Binding_PI3.txt", 
                        SOS1 = "./weights_Binding_SOS.txt", 
                        K55 = "./weights_Binding_K55.txt", 
                        K27 = "./weights_Binding_K27.txt", 
                        K13 = "./weights_Binding_K13.txt", 
                        K19 = "./weights_Binding_K19.txt")
    cat("\n========== Calculating correlation matrix for heatmap ==========\n")
    cor_results <- multimodalallostery_calculate_correlation_matrix(
      input_files = input_files,
      binder_order = BINDER_ORDER,
      binding_sites_map = binding_sites_map
    )
    cor_matrix <- cor_results$cor_matrix
    p_matrix <- cor_results$p_matrix
    cat("\n========== Correlation Matrix ==========\n")
    print(cor_matrix)
    cat("\n========== Generating heatmap with fixed order ==========\n")
    output_heatmap <- "./results/s5a_correlation_heatmap_8binder.pdf"
    multimodalallostery_plot_correlation_heatmap(cor_matrix = cor_matrix, p_matrix = p_matrix, output_file = output_heatmap, width = 10, height = 8)
    cat("\n========== Analysis Complete ==========\n")
    cat("Heatmap saved to:", output_heatmap, "\n")
  })
  invisible(result)
}


run_figure_s5_a()

# ✅
# ---- FigureS5B ----
run_figure_s5_b <- function() {
  result <- local({
    input_file <- "./mutation_classification_summary_8binder.csv"
    output_file <- "./results/s5b_allosteric_only_barplot in pairwise comparison between BI1 and BI2 2.pdf"
    pair_order <- c("RAF1 vs K13", "RAF1 vs K19", "RALGDS vs K13", "RALGDS vs K19", "PI3KCG vs K13", "PI3KCG vs K19", "SOS1 vs K13", 
                    "SOS1 vs K19", "K55 vs K13", "K55 vs K19", "K27 vs K13", "K27 vs K19")
    colors <- c(`Allosteric only in X` = "#1B38A6", `Allosteric only in Y` = "#75C2F6")
    plot_data <- multimodalallostery_prepare_allosteric_pair_counts(data.table::fread(input_file), pair_order)
    p <- multimodalallostery_plot_allosteric_pair_counts(plot_data, colors)
    print(p)
    ggplot2::ggsave(output_file, plot = p, width = 5, height = 3, dpi = 300)
  })
  invisible(result)
}


run_figure_s5_b()


# ✅
# ---- FigureS5C ----
run_figure_s5_c <- function() {
  result <- local({
    binding_sites_map <- list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K27 = c(21, 24, 25, 27, 31, 33, 36, 
                                                                                                    38, 39, 40, 41, 43, 52, 54, 67, 70, 71), K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 
                                                                                                                                                     107, 129, 133, 136, 137, 138), K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 
                                                                                                                                                                                            133, 136, 137))
    NBP_RESIDUES <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    cat("\n========================================================")
    cat("\nSTEP 1 — RAF1 vs K27 anticorrelated mutations")
    cat("\n========================================================\n")
    raf1_file <- "./weights_Binding_RAF.txt"
    k27_file <- "./weights_Binding_K27.txt"
    k13_file <- "./weights_Binding_K13.txt"
    k19_file <- "./weights_Binding_K19.txt"
    anno_file <- "./anno_final_for_8.csv"
    raf1_k27_anticorrelated <- multimodalallostery_get_raf1_vs_k27_anticorrelated(raf1_file, k27_file, anno_file, binding_sites_map, NBP_RESIDUES)
    cat("\n========================================================")
    cat("\nSTEP 2 — Load K13 vs K19 data")
    cat("\n========================================================\n")
    k13_k19_data <- multimodalallostery_load_k13_vs_k19_data(k13_file, k19_file, anno_file, binding_sites_map)
    cat("\n========================================================")
    cat("\nSTEP 3 — Generate plot")
    cat("\n========================================================\n")
    mapping_plot <- multimodalallostery_plot_mapped_anticorrelated(k13_k19_data, raf1_k27_anticorrelated, NBP_RESIDUES, xlim = c(-1.5, 3), ylim = c(-1.5, 3))
    print(mapping_plot)
    output_plot_path <- "./results/s5c_RAF1_K27_anticorrelated_mapped_to_K13_vs_K19.pdf"
    ggplot2::ggsave(output_plot_path, plot = mapping_plot, width = 4.5, height = 5, device = grDevices::cairo_pdf)
    cat("\nPlot saved to:\n")
    cat(output_plot_path, "\n")
  })
  invisible(result)
}

run_figure_s5_c()


# ✅
# ---- FigureS5D_K13 ----
run_figure_s5_d_k13 <- function() {
  result <- local({
    data_dir <- "./"
    tables <- list(RAF1_VS_K27 = data.table::fread(file.path(data_dir, "RAF1_vs_K27_anticorrelated_mutations.csv")), RAF1_VS_K13 = data.table::fread(file.path(data_dir, 
                                                                                                                                                               "RAF1_vs_K13_anticorrelated_mutations.csv")), K27_VS_K13 = data.table::fread(file.path(data_dir, "K27_vs_K13_anticorrelated_mutations.csv")))
    venn_plot <- multimodalallostery_plot_mutation_venn(multimodalallostery_prepare_mutation_sets(tables))
    print(venn_plot)
    ggplot2::ggsave("./results/s5d_RAF1 vs K27_K27 vs K13_RAF1 vs K13_anticorrelated muts_overlap_vennplot.pdf", plot = venn_plot, 
                    width = 5.5, height = 5, dpi = 300, units = "in", device = "pdf")
  })
  invisible(result)
}

run_figure_s5_d_k13()


# 
# ---- FigureS5D_K19 ----
run_figure_s5_d_k19 <- function() {
  result <- local({
    data_dir <- "./"
    tables <- list(RAF1_VS_K27 = data.table::fread(file.path(data_dir, "RAF1_vs_K27_anticorrelated_mutations.csv")), RAF1_VS_K19 = data.table::fread(file.path(data_dir, 
                                                                                                                                                               "RAF1_vs_K19_anticorrelated_mutations.csv")), K27_VS_K19 = data.table::fread(file.path(data_dir, "K27_vs_K19_anticorrelated_mutations.csv")))
    venn_plot <- multimodalallostery_plot_mutation_venn(multimodalallostery_prepare_mutation_sets(tables))
    print(venn_plot)
    ggplot2::ggsave("./results/s5d_RAF1 vs K27_K27 vs K19_RAF1 vs K19_anticorrelated muts_overlap_vennplot.pdf", plot = venn_plot, 
                    width = 5.5, height = 5, dpi = 300, units = "in", device = "pdf")
  })
  invisible(result)
}

run_figure_s5_d_k19()


# ✅
# ---- FigureS5E ----
run_figure_s5_e <- function() {
  result <- local({
    binding_sites_map <- list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K27 = c(21, 24, 25, 27, 31, 33, 36, 
                                                                                                    38, 39, 40, 41, 43, 52, 54, 67, 70, 71), K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 
                                                                                                                                                     107, 129, 133, 136, 137, 138), K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 
                                                                                                                                                                                            133, 136, 137))
    NBP_RESIDUES <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    cat("\n========================================================")
    cat("\nSTEP 1 — Get K13 vs K19 both promoting mutations")
    cat("\n========================================================\n")
    k13_file <- "./weights_Binding_K13.txt"
    k19_file <- "./weights_Binding_K19.txt"
    k27_file <- "./weights_Binding_K27.txt"
    anno_file <- "./anno_final_for_8.csv"
    both_promoting_muts <- multimodalallostery_get_k13_vs_k19_both_promoting(k13_file, k19_file, anno_file, binding_sites_map, NBP_RESIDUES)
    cat("\n========================================================")
    cat("\nSTEP 2 — Load K27 vs K13 data")
    cat("\n========================================================\n")
    k27_k13_data <- multimodalallostery_load_k27_vs_k13_data(k27_file, k13_file, anno_file, binding_sites_map)
    cat("\n========================================================")
    cat("\nSTEP 3 — Load K27 vs K19 data")
    cat("\n========================================================\n")
    k27_k19_data <- multimodalallostery_load_k27_vs_k19_data(k27_file, k19_file, anno_file, binding_sites_map)
    cat("\n========================================================")
    cat("\nSTEP 4 — Generate plot: K27 vs K13 (fixed axis)")
    cat("\n========================================================\n")
    mapping_plot_K27_K13 <- multimodalallostery_plot_mapped_both_promoting_k27_k13(k27_k13_data, both_promoting_muts, xlim = c(-1.5, 3), ylim = c(-1.5, 
                                                                                                                              3))
    print(mapping_plot_K27_K13)
    cat("\n========================================================")
    cat("\nSTEP 5 — Generate plot: K27 vs K19 (fixed axis)")
    cat("\n========================================================\n")
    mapping_plot_K27_K19 <- multimodalallostery_plot_mapped_both_promoting_k27_k19(k27_k19_data, both_promoting_muts, xlim = c(-1.5, 3), ylim = c(-1.5, 
                                                                                                                              3))
    print(mapping_plot_K27_K19)
    output_dir <- "./results/"
    ggplot2::ggsave(file.path(output_dir, "s5e_K13_K19_both_promoting_mapped_to_K27_vs_K13.pdf"), plot = mapping_plot_K27_K13, 
                    width = 4.5, height = 5, device = grDevices::cairo_pdf)
    ggplot2::ggsave(file.path(output_dir, "s5e_K13_K19_both_promoting_mapped_to_K27_vs_K19.pdf"), plot = mapping_plot_K27_K19, 
                    width = 4.5, height = 5, device = grDevices::cairo_pdf)
    cat("\nPlots saved to:\n")
    cat(file.path(output_dir, "s5e_K13_K19_both_promoting_mapped_to_K27_vs_K13.pdf"), "\n")
    cat(file.path(output_dir, "s5e_K13_K19_both_promoting_mapped_to_K27_vs_K19.pdf"), "\n")
  })
  invisible(result)
}

run_figure_s5_e()

# ✅
# ---- FigureS6A ----
run_figure_s6_a <- function() {
  result <- local({
    PAIR_ORDER <- c("RAF1 vs K55")
    binding_sites_map <- list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K55 = c(5, 24, 25, 31, 33, 36, 37, 
                                                                                                    38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74))
    core_residues <- c(4, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 40, 42, 44, 46, 51, 52, 53, 54, 55, 
                       56, 57, 58, 68, 71, 72, 75, 77, 78, 79, 80, 81, 82, 83, 84, 89, 90, 92, 93, 96, 97, 99, 100, 101, 103, 109, 110, 111, 
                       112, 113, 114, 115, 116, 118, 125, 130, 133, 134, 137, 139, 141, 142, 143, 144, 145, 146, 151, 152, 155, 156, 157, 158, 
                       159, 160, 162, 163)
    NBP <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    functional_loop_residues <- c(10, 11, 12, 13, 14, 15, 16, 17, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
                                  40, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76)
    beta_sheets <- c(3:9, 38:44, 51:57, 77:84, 109:115, 139:143)
    surface_residues <- c(1, 2, 3, 5, 12, 13, 25:39, 41, 43, 45, 47:50, 59:67, 69, 70, 73, 74, 76, 85:88, 91, 94, 95, 98, 102, 
                          104:108, 117, 119:124, 126:129, 131, 132, 135, 136, 138, 140, 147:150, 153, 154, 161, 164:188)
    alpha_helices <- c(15:24, 67:73, 87:104, 127:136, 148:166)
    input_files <- list(RAF1 = "./weights_Binding_RAF.txt", 
                        K55 = "./weights_Binding_K55.txt")
    anno <- "./anno_final_for_8.csv"
    structure_regions <- list(Core = core_residues, Surface = surface_residues, NBP = NBP, `Functional Loop` = functional_loop_residues, 
                              `Beta Sheets` = beta_sheets, `Alpha Helices` = alpha_helices)
    target_pair <- "RAF1 vs K55"
    res <- multimodalallostery_run_pair_region_enrichment(target_pair, input_files, anno, structure_regions, binding_sites_map)
    prepared_plot <- multimodalallostery_prepare_region_enrichment_plot_data(res, target_pair, names(structure_regions))
    color_map <- c(Correlated = "#F4AD0C", `Anti-correlated` = "#1B38A6", Other = "grey80")
    p <- multimodalallostery_plot_region_enrichment(prepared_plot, color_map)
    print(p)
    output_dir <- "./results/"
    ggplot2::ggsave(paste0(output_dir, "s6a_barplot_RAF1_vs_K55_core_surface_nbp_fl_bsheet_ahelix_enrichment.pdf"), plot = p, width = 10, 
                    height = 6, device = grDevices::cairo_pdf)
    outputs <- multimodalallostery_format_region_enrichment_outputs(res)
    data.table::fwrite(outputs$odds, paste0(output_dir, "s6a_OR_values_RAF1_vs_K55_all_six_regions.csv"))
    data.table::fwrite(outputs$fractions, paste0(output_dir, "s6a_fraction_summary_RAF1_vs_K55_all_six_regions.csv"))
  })
  invisible(result)
}


run_figure_s6_a()



# ✅
# ---- FigureS6B ----
run_figure_s6_b <- function() {
  result <- local({
    PAIR_ORDER <- c("RAF1 vs K27")
    binding_sites_map <- list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 
                                                                                                                                                         71))
    core_residues <- c(4, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 40, 42, 44, 46, 51, 52, 53, 54, 55, 
                       56, 57, 58, 68, 71, 72, 75, 77, 78, 79, 80, 81, 82, 83, 84, 89, 90, 92, 93, 96, 97, 99, 100, 101, 103, 109, 110, 111, 
                       112, 113, 114, 115, 116, 118, 125, 130, 133, 134, 137, 139, 141, 142, 143, 144, 145, 146, 151, 152, 155, 156, 157, 158, 
                       159, 160, 162, 163)
    NBP <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    functional_loop_residues <- c(10, 11, 12, 13, 14, 15, 16, 17, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
                                  40, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76)
    beta_sheets <- c(3:9, 38:44, 51:57, 77:84, 109:115, 139:143)
    surface_residues <- c(1, 2, 3, 5, 12, 13, 25:39, 41, 43, 45, 47:50, 59:67, 69, 70, 73, 74, 76, 85:88, 91, 94, 95, 98, 102, 
                          104:108, 117, 119:124, 126:129, 131, 132, 135, 136, 138, 140, 147:150, 153, 154, 161, 164:188)
    alpha_helices <- c(15:24, 67:73, 87:104, 127:136, 148:166)
    input_files <- list(RAF1 = "./weights_Binding_RAF.txt", 
                        K27 = "./weights_Binding_K27.txt")
    anno <- "./anno_final_for_8.csv"
    structure_regions <- list(Core = core_residues, Surface = surface_residues, NBP = NBP, `Functional Loop` = functional_loop_residues, 
                              `Beta Sheets` = beta_sheets, `Alpha Helices` = alpha_helices)
    target_pair <- "RAF1 vs K27"
    res <- multimodalallostery_run_pair_region_enrichment(target_pair, input_files, anno, structure_regions, binding_sites_map)
    prepared_plot <- multimodalallostery_prepare_region_enrichment_plot_data(res, target_pair, names(structure_regions))
    color_map <- c(Correlated = "#F4AD0C", `Anti-correlated` = "#1B38A6", Other = "grey80")
    p <- multimodalallostery_plot_region_enrichment(prepared_plot, color_map)
    print(p)
    output_dir <- "./results/"
    ggplot2::ggsave(paste0(output_dir, "s6b_barplot_RAF1_vs_K27_core_surface_nbp_fl_bsheet_ahelix_enrichment.pdf"), plot = p, width = 10, 
                    height = 6, device = grDevices::cairo_pdf)
    outputs <- multimodalallostery_format_region_enrichment_outputs(res)
    data.table::fwrite(outputs$odds, paste0(output_dir, "s6b_OR_values_RAF1_vs_K27_all_six_regions.csv"))
    data.table::fwrite(outputs$fractions, paste0(output_dir, "s6b_fraction_summary_RAF1_vs_K27_all_six_regions.csv"))
  })
  invisible(result)
}


run_figure_s6_b()


# 
# ---- FigureS6C ----
run_figure_s6_c <- function() {
  result <- local({
    PAIR_ORDER <- c("K27 vs K13")
    binding_sites_map <- list(K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71), K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138))
    core_residues <- c(4, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 40, 42, 44, 46, 51, 52, 53, 54, 55, 
                       56, 57, 58, 68, 71, 72, 75, 77, 78, 79, 80, 81, 82, 83, 84, 89, 90, 92, 93, 96, 97, 99, 100, 101, 103, 109, 110, 111, 
                       112, 113, 114, 115, 116, 118, 125, 130, 133, 134, 137, 139, 141, 142, 143, 144, 145, 146, 151, 152, 155, 156, 157, 158, 
                       159, 160, 162, 163)
    NBP <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    functional_loop_residues <- c(10, 11, 12, 13, 14, 15, 16, 17, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
                                  40, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76)
    beta_sheets <- c(3:9, 38:44, 51:57, 77:84, 109:115, 139:143)
    surface_residues <- c(1, 2, 3, 5, 12, 13, 25:39, 41, 43, 45, 47:50, 59:67, 69, 70, 73, 74, 76, 85:88, 91, 94, 95, 98, 102, 
                          104:108, 117, 119:124, 126:129, 131, 132, 135, 136, 138, 140, 147:150, 153, 154, 161, 164:188)
    alpha_helices <- c(15:24, 67:73, 87:104, 127:136, 148:166)
    input_files <- list(K27 = "./weights_Binding_K27.txt", 
                        K13 = "./weights_Binding_K13.txt")
    anno <- "./anno_final_for_8.csv"
    structure_regions <- list(Core = core_residues, Surface = surface_residues, NBP = NBP, `Functional Loop` = functional_loop_residues, 
                              `Beta Sheets` = beta_sheets, `Alpha Helices` = alpha_helices)
    target_pair <- "K27 vs K13"
    res <- multimodalallostery_run_pair_region_enrichment(target_pair, input_files, anno, structure_regions, binding_sites_map)
    prepared_plot <- multimodalallostery_prepare_region_enrichment_plot_data(res, target_pair, names(structure_regions))
    color_map <- c(Correlated = "#F4AD0C", `Anti-correlated` = "#1B38A6", Other = "grey80")
    p <- multimodalallostery_plot_region_enrichment(prepared_plot, color_map)
    print(p)
    output_dir <- "./results/"
    ggplot2::ggsave(paste0(output_dir, "s6c_barplot_K27_vs_K13_core_surface_nbp_fl_bsheet_ahelix_enrichment.pdf"), plot = p, width = 10, height = 6, 
                    device = grDevices::cairo_pdf)
    outputs <- multimodalallostery_format_region_enrichment_outputs(res)
    data.table::fwrite(outputs$odds, paste0(output_dir, "s6c_OR_values_K27_vs_K13_all_six_regions.csv"))
    data.table::fwrite(outputs$fractions, paste0(output_dir, "s6c_fraction_summary_K27_vs_K13_all_six_regions.csv"))
  })
  invisible(result)
}


run_figure_s6_c()



# ✅
# ---- FigureS6D ----
run_figure_s6_d <- function() {
  result <- local({
    PAIR_ORDER <- c("K13 vs K19")
    binding_sites_map <- list(K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138), K19 = c(68,87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137))
    core_residues <- c(4, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 40, 42, 44, 46, 51, 52, 53, 54, 55, 
                       56, 57, 58, 68, 71, 72, 75, 77, 78, 79, 80, 81, 82, 83, 84, 89, 90, 92, 93, 96, 97, 99, 100, 101, 103, 109, 110, 111, 
                       112, 113, 114, 115, 116, 118, 125, 130, 133, 134, 137, 139, 141, 142, 143, 144, 145, 146, 151, 152, 155, 156, 157, 158, 
                       159, 160, 162, 163)
    NBP <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    functional_loop_residues <- c(10, 11, 12, 13, 14, 15, 16, 17, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
                                  40, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76)
    beta_sheets <- c(3:9, 38:44, 51:57, 77:84, 109:115, 139:143)
    surface_residues <- c(1, 2, 3, 5, 12, 13, 25:39, 41, 43, 45, 47:50, 59:67, 69, 70, 73, 74, 76, 85:88, 91, 94, 95, 98, 102, 
                          104:108, 117, 119:124, 126:129, 131, 132, 135, 136, 138, 140, 147:150, 153, 154, 161, 164:188)
    alpha_helices <- c(15:24, 67:73, 87:104, 127:136, 148:166)
    input_files <- list(K13 = "./weights_Binding_K13.txt", 
                        K19 = "./weights_Binding_K19.txt")
    anno <- "./anno_final_for_8.csv"
    structure_regions <- list(Core = core_residues, Surface = surface_residues, NBP = NBP, `Functional Loop` = functional_loop_residues, 
                              `Beta Sheets` = beta_sheets, `Alpha Helices` = alpha_helices)
    target_pair <- "K13 vs K19"
    res <- multimodalallostery_run_pair_region_enrichment(target_pair, input_files, anno, structure_regions, binding_sites_map)
    prepared_plot <- multimodalallostery_prepare_region_enrichment_plot_data(res, target_pair, names(structure_regions))
    color_map <- c(Correlated = "#F4AD0C", `Anti-correlated` = "#1B38A6", Other = "grey80")
    p <- multimodalallostery_plot_region_enrichment(prepared_plot, color_map)
    print(p)
    output_dir <- "./results/"
    ggplot2::ggsave(paste0(output_dir, "s6d_barplot_K13 VS K19_core_surface_nbp_fl_bsheet_ahelix_enrichment.pdf"), plot = p, width = 10, height = 6, 
                    device = grDevices::cairo_pdf)
    outputs <- multimodalallostery_format_region_enrichment_outputs(res)
    data.table::fwrite(outputs$odds, paste0(output_dir, "s6d_OR_values_K13_vs_K19_all_six_regions.csv"))
    data.table::fwrite(outputs$fractions, paste0(output_dir, "s6d_fraction_summary_K13_vs_K19_all_six_regions.csv"))
  })
  invisible(result)
}


run_figure_s6_d()


# ✅
# ---- FigureS6E ----
run_figure_s6_e <- function() {
  result <- local({
    PAIR_ORDER <- c("RAF1 vs K13")
    binding_sites_map <- list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138))
    input_files <- list(RAF1 = "./weights_Binding_RAF.txt", 
                        K13 = "./weights_Binding_K13.txt")
    anno <- "./anno_final_for_8.csv"
    core_residues <- c(4, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 40, 42, 44, 46, 51, 52, 53, 54, 55, 
                       56, 57, 58, 68, 71, 72, 75, 77, 78, 79, 80, 81, 82, 83, 84, 89, 90, 92, 93, 96, 97, 99, 100, 101, 103, 109, 110, 111, 
                       112, 113, 114, 115, 116, 118, 125, 130, 133, 134, 137, 139, 141, 142, 143, 144, 145, 146, 151, 152, 155, 156, 157, 158, 
                       159, 160, 162, 163)
    NBP <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    functional_loop_residues <- c(10, 11, 12, 13, 14, 15, 16, 17, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 
                                  40, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76)
    beta_sheets <- c(3:9, 38:44, 51:57, 77:84, 109:115, 139:143)
    surface_residues <- c(1, 2, 3, 5, 12, 13, 25:39, 41, 43, 45, 47:50, 59:67, 69, 70, 73, 74, 76, 85:88, 91, 94, 95, 98, 102, 
                          104:108, 117, 119:124, 126:129, 131, 132, 135, 136, 138, 140, 147:150, 153, 154, 161, 164:188)
    alpha_helices <- c(15:24, 67:73, 87:104, 127:136, 148:166)
    structure_regions <- list(Core = core_residues, Surface = surface_residues, NBP = NBP, `Functional Loop` = functional_loop_residues, 
                              `Beta Sheets` = beta_sheets, `Alpha Helices` = alpha_helices)
    contact_shell_file <- "./5binder_contact_shell2.csv"
    contact_shell <- data.table::fread(contact_shell_file)
    second_shell_map <- multimodalallostery_prepare_second_shell_map(contact_shell, names(input_files))
    target_pair <- "RAF1 vs K13"
    res <- multimodalallostery_run_pair_second_shell_enrichment(target_pair, input_files, anno, structure_regions, second_shell_map, binding_sites_map)
    prepared_plot <- multimodalallostery_prepare_allosteric_region_plot_data(res, target_pair, names(structure_regions))
    color_map <- c(`Allosteric only in X` = "#C68EFD", `Allosteric only in Y` = "#09B636", Other = "grey80")
    p <- multimodalallostery_plot_allosteric_region_enrichment(prepared_plot, color_map)
    print(p)
    output_dir <- "./results/"
    ggplot2::ggsave(paste0(output_dir, "s6e_barplot_RAF1_vs_K13_allosteric_only_with_second_shell2.pdf"), plot = p, width = 10, height = 6, 
                    device = grDevices::cairo_pdf)
    or_output <- multimodalallostery_format_allosteric_region_outputs(res)
    data.table::fwrite(or_output, paste0(output_dir, "s6e_OR_values_RAF1_vs_K13_allosteric_only_with_second_shell.csv"))
  })
  invisible(result)
}

run_figure_s6_e()


# ✅
# ---- FigureS6F ----
run_figure_s6_f <- function() {
  result <- local({
    input_file <- "./RAF1_vs_K13_correlation_results_with_WT_0.csv"
    output_file <- "./results/s6f_PanelC_PearsonR_annotation5.pdf"
    annotation_tracks <- list(BI1 = c(5, 21, 24, 25, 27, 29, 31, 33, 36, 37, 38, 39, 40, 41, 43, 52, 54, 56, 64, 66, 67, 70, 
                                      71, 73, 74), BI2 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 108, 125, 129, 133, 
                                                           136, 137, 138), `β-strand 1` = 3:9, `β-strand 2` = 38:44, `β-strand 3` = 51:57, `β-strand 4` = 77:84, `β-strand 5` = 109:115, 
                              `β-strand 6` = 139:143, NBP = c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 
                                                              145, 146, 147), `γ-phosphate contact` = c(16, 35, 60, 12, 13, 34))
    plot_df <- multimodalallostery_prepare_residue_correlations(data.table::fread(input_file))
    annotation_df <- multimodalallostery_prepare_annotation_tracks(annotation_tracks, plot_df)
    panel <- patchwork::wrap_plots(multimodalallostery_plot_residue_correlations(plot_df), multimodalallostery_plot_annotation_tracks(annotation_df, levels(plot_df$residue_order)), 
                                   ncol = 1, heights = c(4, 2.5))
    panel
    ggplot2::ggsave(output_file, panel, width = 16, height = 6, dpi = 600)
  })
  invisible(result)
}

run_figure_s6_f()


# ✅
# ---- FigureS6G ----
run_figure_s6_g <- function() {
  result <- local({
    COLOR_MAP <- c(`Not significant (FDR >= 0.05)` = "grey90", `Other (neutral in both)` = "grey90", `Both promoting` = "#FFB0A5", 
                   `Both disrupting` = "#F4270C", `Promoting in X / Disrupting in Y` = "#F4AD0C", `Disrupting in X / Promoting in Y` = "#F1DD10", 
                   `Allosteric only in X` = "#1B38A6", `Allosteric only in Y` = "#75C2F6")
    LEGEND_ORDER <- c("Both promoting", "Both disrupting", "Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y", 
                      "Allosteric only in X", "Allosteric only in Y", "Other (neutral in both)", "Not significant (FDR >= 0.05)")
    binding_sites_map <- list(RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71), K55 = c(5, 24, 25, 31, 33, 36, 37, 
                                                                                                    38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74), K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 
                                                                                                                                                         71), K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138), K19 = c(68, 
                                                                                                                                                                                                                                                                                     87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137))
    NBP_residues <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
    Switch_I_residues <- c(25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40)
    Switch_II_residues <- c(58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76)
    anno <- "./anno_final_for_8.csv"
    input_files <- list(RAF1 = "./weights_Binding_RAF.txt", 
                        K13 = "./weights_Binding_K13.txt")
    full_result <- multimodalallostery_analyze_all_sites(input_x = input_files$RAF1, input_y = input_files$K13, assay_x = "RAF1", assay_y = "K13", 
                                     anno = anno, legend_order = LEGEND_ORDER, nbp_residues = NBP_residues, switch_i_residues = Switch_I_residues, switch_ii_residues = Switch_II_residues)
    cat("\n=== Data Validation ===\n")
    cat("Total number of rows:", nrow(full_result$data), "\n")
    cat("WT row number:", nrow(full_result$data[mt_codon == "WT"]), "\n")
    cat("Number of sites:", length(unique(full_result$data$Pos_real)), "\n")
    p_pos145 <- multimodalallostery_plot_site_scatter(full_analysis_result = full_result, target_position = 145, point_size = 5, base_size = 20, 
                                  show_labels = TRUE, label_all = TRUE, xlim = c(-1.5, 3.3), ylim = c(-1.5, 3.3), show_WT = TRUE, WT_point_size = 5, WT_color = "black", 
                                  WT_shape = 19, WT_label = "WT")
    ggplot2::ggsave("./results/RAF1_vs_K13_scatter_plot_position_145_with_WT.pdf", 
                    plot = p_pos145, width = 6.5, height = 8, device = grDevices::cairo_pdf)
    p15 <- multimodalallostery_plot_site_scatter(full_analysis_result = full_result, target_position = 15, point_size = 5, base_size = 20, show_labels = TRUE, 
                             label_all = TRUE, xlim = c(-1.5, 3.3), ylim = c(-1.5, 3.3), show_WT = TRUE, WT_point_size = 5, WT_color = "black", WT_shape = 19, 
                             WT_label = "WT")
    ggplot2::ggsave("./results/RAF1_vs_K13_scatter_plot_position_15_with_WT.pdf", 
                    plot = p15, width = 6.5, height = 8, device = grDevices::cairo_pdf)
    p48 <- multimodalallostery_plot_site_scatter(full_analysis_result = full_result, target_position = 48, point_size = 5, base_size = 20, show_labels = TRUE, 
                             label_all = TRUE, xlim = c(-1.5, 3.3), ylim = c(-1.5, 3.3), show_WT = TRUE, WT_point_size = 5, WT_color = "black", WT_shape = 19, 
                             WT_label = "WT")
    ggplot2::ggsave("./results/RAF1_vs_K13_scatter_plot_position_48_with_WT.pdf", 
                    plot = p48, width = 6.5, height = 8, device = grDevices::cairo_pdf)
    p55 <- multimodalallostery_plot_site_scatter(full_analysis_result = full_result, target_position = 55, point_size = 5, base_size = 20, show_labels = TRUE, 
                             label_all = TRUE, xlim = c(-1.5, 3.3), ylim = c(-1.5, 3.3), show_WT = TRUE, WT_point_size = 5, WT_color = "black", WT_shape = 19, 
                             WT_label = "WT")
    ggplot2::ggsave("./results/RAF1_vs_K13_scatter_plot_position_55_with_WT.pdf", 
                    plot = p55, width = 6.5, height = 8, device = grDevices::cairo_pdf)
    p77 <- multimodalallostery_plot_site_scatter(full_analysis_result = full_result, target_position = 77, point_size = 5, base_size = 20, show_labels = TRUE, 
                             label_all = TRUE, xlim = c(-1.5, 3.3), ylim = c(-1.5, 3.3), show_WT = TRUE, WT_point_size = 5, WT_color = "black", WT_shape = 19, 
                             WT_label = "WT")
    ggplot2::ggsave("./results/RAF1_vs_K13_scatter_plot_position_77_with_WT.pdf", 
                    plot = p77, width = 6.5, height = 8, device = grDevices::cairo_pdf)
    p163 <- multimodalallostery_plot_site_scatter(full_analysis_result = full_result, target_position = 163, point_size = 5, base_size = 20, show_labels = TRUE, 
                              label_all = TRUE, xlim = c(-1.5, 3.3), ylim = c(-1.5, 3.3), show_WT = TRUE, WT_point_size = 5, WT_color = "black", WT_shape = 19, 
                              WT_label = "WT")
    ggplot2::ggsave("./results/RAF1_vs_K13_scatter_plot_position_163_with_WT.pdf", 
                    plot = p163, width = 6.5, height = 8, device = grDevices::cairo_pdf)
  })
  invisible(result)
}


run_figure_s6_g()


run_all_figures <- function(panels = NULL) {
  runners <- list(
    Figure1E = run_figure1_e,
    Figure1F = run_figure1_f,
    Figure1G = run_figure1_g,
    Figure1H = run_figure1_h,
    Figure1I = run_figure1_i,
    Figure2A = run_figure2_a,
    Figure2C = run_figure2_c,
    Figure2D = run_figure2_d,
    Figure3A = run_figure3_a,
    Figure3C = run_figure3_c,
    Figure3D = run_figure3_d,
    Figure3E = run_figure3_e,
    Figure3F = run_figure3_f,
    Figure4B = run_figure4_b,
    Figure5B = run_figure5_b,
    Figure6B = run_figure6_b,
    FigureS1B = run_figure_s1_b,
    FigureS1C = run_figure_s1_c,
    FigureS1D = run_figure_s1_d,
    FigureS2A = run_figure_s2_a,
    FigureS2B = run_figure_s2_b,
    FigureS2C = run_figure_s2_c,
    FigureS2D = run_figure_s2_d,
    FigureS2F = run_figure_s2_f,
    FigureS3A = run_figure_s3_a,
    FigureS3B = run_figure_s3_b,
    FigureS3C = run_figure_s3_c,
    FigureS3E = run_figure_s3_e,
    FigureS3F = run_figure_s3_f,
    FigureS3G = run_figure_s3_g,
    FigureS3H = run_figure_s3_h,
    FigureS4A = run_figure_s4_a,
    FigureS5A = run_figure_s5_a,
    FigureS5B = run_figure_s5_b,
    FigureS5C = run_figure_s5_c,
    FigureS5D_K13 = run_figure_s5_d_k13,
    FigureS5D_K19 = run_figure_s5_d_k19,
    FigureS5E = run_figure_s5_e,
    FigureS6A = run_figure_s6_a,
    FigureS6B = run_figure_s6_b,
    FigureS6C = run_figure_s6_c,
    FigureS6D = run_figure_s6_d,
    FigureS6E = run_figure_s6_e,
    FigureS6F = run_figure_s6_f,
    FigureS6G = run_figure_s6_g
  )
  if (!is.null(panels)) runners <- runners[panels]
  lapply(runners, function(run_panel) run_panel())
}
