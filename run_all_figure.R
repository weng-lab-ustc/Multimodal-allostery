# Formal September 2026 figure workflow for multimodalallostery.
#
# This orchestration script follows the explicit Figure-level layout of
# package_refactor_20260904/run_all_figure.R. Scientific calculations remain in
# the package; paths, panel order, package calls, and outputs are defined here.

library(krasddpcams)
library(data.table)
library(wlab.block)

#options(timeout = 600)
#remotes::install_github(
#  "weng-lab-ustc/Multimodal-allostery",
#  dependencies = TRUE
#)


library(multimodalallostery)

setwd("~/Library/CloudStorage/OneDrive-个人/文档/Script/multimodalallostery/")


# ---- Configuration: edit paths here when data move --------------------------
WEIGHTS_DIR <- "./Files for plot of multimodalallostery/Energy data/"
PREDICTIONS_FILE <- "./Files for plot of multimodalallostery/Energy data/predicted_phenotypes_all.txt"
ANNOTATION_5 <- "./Files for plot of multimodalallostery/anno_final_for_5.csv"
ANNOTATION_8 <- "./Files for plot of multimodalallostery/anno_final_for_8.csv"
CONTACT_SHELL <- "./Files for plot of multimodalallostery/5binder_contact_shell2.csv"
FITNESS_DIR <- "./Files for plot of multimodalallostery/fitness_RData/"
FITNESS_MERGE_DIR <- "./Files for plot of multimodalallostery/fitness_RData_merge_version/"
MANUSCRIPT_FIGURE_DIR <- "./Files for plot of multimodalallostery/"
OUTPUT_DIR <- "./Files for plot of multimodalallostery/multimodalallostery_results2"


# Set to NULL to run all 33 panels in order, or to a subset such as
# c("F1e", "F1f"). Sourcing this file starts the selected workflow.
PANELS_TO_RUN <- NULL
OVERWRITE_OUTPUTS <- FALSE

# Ensure Chinese path names work when Rscript inherits a non-Windows locale.
invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))

message("Weights input: ", WEIGHTS_DIR)
message("Fitness input: ", FITNESS_DIR)
message("Fitness merge input: ", FITNESS_MERGE_DIR)
message("Annotation input: ", ANNOTATION_8)
message("Output directory: ", OUTPUT_DIR)

steps <- data.frame(
  panel = c("F1e", "F1f", "F1g", "F1h", "F1i", "F2a", "F2c", "F2d",
            "F3a", "F3c", "F3d", "F3e", "F3f", "F4b", "F5a", "F6b",
            "SF1b", "SF1c", "SF1d",
            "SF2a", "SF2b", "SF2c", "SF2d", "SF2e", "SF2f", "SF2g",
            "SF3a", "SF3b", "SF3c", "SF3e", "SF3f", "SF3g", "SF3h",
            "SF4a", "SF5a", "SF5b", "SF5c", "SF5d", "SF5e",
            "SF6a", "SF6b", "SF6c", "SF6d"),
  group = c(rep("Figure1", 5), rep("Figure2", 3), rep("Figure3", 5),
            "Figure4", "Figure5", "Figure6", rep("SupplementaryFigure1", 3),
            rep("SupplementaryFigure2", 7),rep("SupplementaryFigure3", 7),
            "SupplementaryFigure4", rep("SupplementaryFigure5", 5),
            rep("SupplementaryFigure6", 4)),
  stringsAsFactors = FALSE
)

weight <- function(assay) {
  file_assay <- if (identical(assay, "RAF1")) "RAF" else assay
  path <- file.path(WEIGHTS_DIR, paste0("weights_Binding_", file_assay, ".txt"))
  require_file(path, paste0("binding weights for ", assay))
}
folding_weight <- function() {
  require_file(file.path(WEIGHTS_DIR, "weights_Folding.txt"), "folding weights")
}
require_file <- function(path, label) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path) || !file.exists(path)) {
    stop("Missing ", label, ": ", path)
  }
  path
}

input <- function(key) {
  paths <- list(annotation_5 = ANNOTATION_5, annotation_8 = ANNOTATION_8,
                contact_shell = CONTACT_SHELL)
  require_file(paths[[key]], key)
}
save_pdf <- function(plot, panel, label, width = 6, height = 5) {
  path <- file.path(OUTPUT_DIR, steps$group[match(panel, steps$panel)],
                    paste0(panel, "_", label, ".pdf"))
  if (file.exists(path) && !OVERWRITE_OUTPUTS) stop("Output exists: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height,
                  device = grDevices::cairo_pdf)
  path
}
save_csv <- function(data, panel, label) {
  path <- file.path(OUTPUT_DIR, steps$group[match(panel, steps$panel)],
                    paste0(panel, "_", label, ".csv"))
  if (file.exists(path) && !OVERWRITE_OUTPUTS) stop("Output exists: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data, path, row.names = FALSE)
  path
}

defs <- multimodalallostery_structure_definitions()

### This dataset is used to plot fitness density (with different types of variants merged into a single data frame).
fitness_merge_blocks <- list(
  Abundance = file.path(FITNESS_MERGE_DIR, c(
    "CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
    "MA_RAS_abundance_2_fitness_replicates_fullseq_merge.RData",
    "MA_RAS_abundance_3_fitness_replicates_fullseq_merge.RData")),
  K13 = file.path(FITNESS_MERGE_DIR, c(
    "MA_RAS_binding_K13_1_fitness_replicates_fullseq_merge.RData",
    "MA_RAS_binding_K13_2_fitness_replicates_fullseq_merge.RData",
    "MA_RAS_binding_K13_3_fitness_replicates_fullseq_merge.RData")),
  K19 = file.path(FITNESS_MERGE_DIR, c(
    "MA_RAS_binding_K19_1_fitness_replicates_fullseq_merge.RData",
    "MA_RAS_binding_K19_2_fitness_replicates_fullseq_merge.RData",
    "MA_RAS_binding_K19_3_fitness_replicates_fullseq_merge.RData"))
)

fitness_blocks <- list(
  Abundance = file.path(FITNESS_DIR, c(
    "CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
    "MA_RAS_abundance_2_fitness_replicates_fullseq.RData",
    "MA_RAS_abundance_3_fitness_replicates_fullseq.RData")),
  K13 = file.path(FITNESS_DIR, c(
    "MA_RAS_binding_K13_1_fitness_replicates_fullseq.RData",
    "MA_RAS_binding_K13_2_fitness_replicates_fullseq.RData",
    "MA_RAS_binding_K13_3_fitness_replicates_fullseq.RData")),
  K19 = file.path(FITNESS_DIR, c(
    "MA_RAS_binding_K19_1_fitness_replicates_fullseq.RData",
    "MA_RAS_binding_K19_2_fitness_replicates_fullseq.RData",
    "MA_RAS_binding_K19_3_fitness_replicates_fullseq.RData")),
  RAF1 = file.path(FITNESS_DIR, c(
    "CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
    "MA_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
    "MA_RAS_binding_RAF_3_fitness_replicates_fullseq.RData"))
)

# The scripts use the KRAS sequence without the initiating methionine for
# mutation-indexed heatmaps and with it for the sequence annotation panel.
kras <- "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
wt_aa <- substring(kras, 2L)
sheet_rects <- data.frame(xstart = c(3, 38, 51, 77, 109, 139),
                          xend = c(9, 44, 57, 84, 115, 143),
                          col = paste0("b", 1:6))
helix_rects <- data.frame(xstart = c(15, 67, 87, 127, 148),
                          xend = c(24, 73, 104, 136, 166),
                          col = paste0("alpha", 1:5))
# Figure 1: fitness, energy, and sequence context.
## F1e
run_F1e <- function() {
  for (assay in names(fitness_merge_blocks)) {
    files <- fitness_merge_blocks[[assay]]
    invisible(lapply(files, require_file, label = paste0(assay, " fitness block")))
    p <- multimodalallostery_plot_fitness_density(
      assay_type = assay, block1 = files[1], block2 = files[2], block3 = files[3])
    save_pdf(p, "F1e", assay)
  }
}
run_F1e()


## F1f
run_F1f <- function() {
  anno <- data.table::fread(input("annotation_5"))
  anno[, Pos_real := Pos] 
  summery <- ddG_data_assay(input = file.path(WEIGHTS_DIR, "weights_Binding_RAF.txt"),wt_aa = wt_aa)
  anno <- merge(anno,summery,by = "Pos_real",all = TRUE)
  abundance_files <- fitness_blocks$Abundance
  invisible(lapply(abundance_files, require_file, label = "Abundance fitness block"))
  stab <- multimodalallostery_normalize_fitness(
    abundance_files[1], abundance_files[2], abundance_files[3])
  K13 <- multimodalallostery_normalize_fitness(
    fitness_blocks$K13[1], fitness_blocks$K13[2], fitness_blocks$K13[3])
  K19 <- multimodalallostery_normalize_fitness(
    fitness_blocks$K19[1], fitness_blocks$K19[2], fitness_blocks$K19[3])
  RAF1 <- multimodalallostery_normalize_fitness(
    fitness_blocks$RAF1[1], fitness_blocks$RAF1[2], fitness_blocks$RAF1[3])
  for (assay in c("K13", "K19", "RAF1")) {
    files <- fitness_blocks[[assay]]
    invisible(lapply(files, require_file, label = paste0(assay, " fitness block")))
    merged <- switch(assay,
                     K13 = multimodalallostery_merge_dimsum_data(stab, K13),
                     K19 = multimodalallostery_merge_dimsum_data(stab, K19),
                     RAF1 = multimodalallostery_merge_dimsum_data(stab, RAF1))
    dat <- multimodalallostery_identify_mutation_positions(merged, wt_aa)
    save_pdf(multimodalallostery_plot_binding_fitness(dat, assay, anno), "F1f", assay)
  }
}

run_F1f()

## F1g
run_F1g <- function() {
  anno <- data.table::fread(input("annotation_5"))
  for (assay in c("K13", "K19", "RAF1")) {
    p <- multimodalallostery_plot_scatter_dd_gb_dd_gf( folding_weight(), "Folding",
                                                       weight(assay), assay, anno, binder = assay,
                                                       colour_scheme = c("#F4270C", "#1B38A6"))
    save_pdf(p, "F1g", assay)
  }
}

run_F1g()

## F1h
run_F1h <- function() {
  assays <- c("Folding", "RAF", "K13", "K19")
  for (assay in assays) {
    path <- if (assay == "Folding") folding_weight() else weight(assay)
    save_pdf(multimodalallostery_dd_g_heatmap( path, wt_aa, title = assay), "F1h", assay,
             width = 20, height = 6)
  }
}

run_F1h()

## F1i
run_F1i <- function() {
  prepared <- multimodalallostery_prepare_sequence_annotation( kras,
                       c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138),
                       c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71))
  save_pdf(multimodalallostery_plot_sequence_annotation( prepared,
                    c(12:18, 28:30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145:147),
                    data.frame(xstart = c(10, 25, 58), xend = c(17, 40, 76),
                               col = c("P-loop", "switch I", "switch II")),
                    defs$core_residues, sheet_rects, helix_rects),
           "F1i", "sequence", width = 18, height = 4)
}

run_F1i()


# Figure 2: interface residue effects.

## F2a
run_F2a <- function() {
  sites <- list(K13 = c(63, 105, 106, 98, 107, 101, 102, 99, 136, 95, 137, 94, 133, 90, 129, 87, 91, 88),
                K19 = c(98, 107, 101, 102, 99, 136, 95, 137, 94, 133, 90, 129, 87, 91, 88, 68, 108))
  for (assay in names(sites)) {
    labels <- stats::setNames(as.character(sites[[assay]]), sites[[assay]])
    p <- multimodalallostery_plot_binding_interface_residue_median_dd_g_heatmap(
                  weight(assay), sites[[assay]], labels, title = assay)
    save_pdf(p, "F2a", assay, width = 12, height = 6)
  }
}

run_F2a()


## F2c
run_F2c <- function() {
  residues <- c("K88", "E91", "T87", "Q129", "F90", "L133", "H94", "Y137",
                "H95", "R68", "S136", "Q99", "R102", "K101", "E107", "E98")
  for (assay in c("K13", "K19")) {
    out <- file.path(OUTPUT_DIR, "Figure2", paste0("F2c_", assay, ".pdf"))
    if (file.exists(out)) stop("Output exists: ", out)
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    multimodalallostery_plot_dd_g_beeswarm( weight(assay), assay, residues, out)
  }
}

run_F2c()


## F2d
run_F2d <- function() {
  a <- multimodalallostery_load_dd_g_data( weight("K13"), "K13")
  b <- multimodalallostery_load_dd_g_data( weight("K19"), "K19")
  dat <- multimodalallostery_prepare_interface_ddg_comparison( a, b,
                  c(98, 107, 101, 102, 99, 136, 95, 137, 94, 133, 90, 129, 87, 91, 88),
                  c(98, 107, 101, 102, 99, 136, 95, 137, 94, 133, 90, 129, 87, 91, 88),
                  "K13", "K19")
  analysis <- multimodalallostery_analyze_correlation_outliers( dat, "K13", "K19", 10)
  save_pdf(multimodalallostery_plot_ddg_correlation_outliers( analysis, "K13", "K19"), "F2d", "interface")
}

run_F2d()


# Figure 3: single-assay allostery and spatial decay.

## F3a
run_F3a <- function() {
  for (assay in c("K13", "RAF1")) {
    result <- multimodalallostery_manhatta_plot_single_assay( weight(assay), assay,
                       input("annotation_5"), sheet_rects, helix_rects, wt_aa,
                       threshold = 0.4)
    save_pdf(result$plot, "F3a", assay, width = 16, height = 6)
  }
}

run_F3a()

## F3c
run_F3c <- function() for (assay in c("RAF1", "K13"))
  save_pdf(multimodalallostery_plot_energy_distance_decay_expfit( weight(assay),
                    assay, input("annotation_5")), "F3c", assay)

run_F3c()

## F3d
run_F3d <- function() for (assay in c("RAF1", "K13"))
  save_pdf(multimodalallostery_plot_energy_distance_decay_expfit_directional( weight(assay),
                    assay, input("annotation_5")), "F3d", assay)
run_F3d()

## F3e
run_F3e <- function() for (assay in c("RAF1", "K13"))
  save_pdf(multimodalallostery_plot_energy_distance_decay_expfit_contact_shell( weight(assay),
                    assay, data.table::fread(input("contact_shell"))), "F3e", assay)
run_F3e()

## F3f
run_F3f <- function() for (assay in c("RAF1", "K13"))
  save_pdf(multimodalallostery_plot_energy_distance_decay_expfit_contact_shell_directional( weight(assay),
                    assay, data.table::fread(input("contact_shell"))), "F3f", assay)
run_F3f()

# Figure 4: shared allosteric hotspots; 

## F4b
run_F4b <- function() {
  anno <- data.table::fread(input("annotation_5"))
  colour_scheme <- list(blue = "#1B38A6", red = "#F4270C")
  allosteric_list <- list(K13 = c(10, 145, 151),K19 = c(10, 145, 151),RAF1 = c(10, 15, 16, 17, 20, 22, 28, 32, 34, 35,54, 55, 57, 58, 59, 60, 77, 144, 145, 146, 163))
  save_pdf(multimodalallostery_plot_triple_dd_g_heatmap(weight("K13"), weight("K19"), weight("RAF"), anno, wt_aa,
      colour_scheme = colour_scheme,allosteric_sites_list = allosteric_list,legend_limits = c(-1.3, 3)),"F4b", "hotspots", 20, 6)
}

run_F4b()

# Figure 5: pairwise mutation classes;
## F5
run_F5 <- function() {
  pairs <- list(
    c("RAF1", "K13"), c("RAF1", "K19"), c("RAF1", "K55"),
    c("RAF1", "K27"), c("K27", "K13"), c("K19", "K13")
  )
  for (pair in pairs) {x <- pair[1]; y <- pair[2]
  result <- multimodalallostery_analyze_protein_pair(x, y, weight(x), weight(y), input("annotation_5"),fixed_threshold = 0.4)
  save_pdf(multimodalallostery_plot_protein_pair(result),"F5a", paste0(x, "_", y), 4.5, 6)
  save_csv(result$data,"F5a",paste0(x, "_", y))}
}

run_F5()


# Figure 6: structural partition of those classes.
## F6b
run_F6b <- function() {
  result <- multimodalallostery_analyze_protein_pair("RAF1","K13",weight("RAF"),weight("K13"),
                                                     input("annotation_8"),fixed_threshold = 0.40)
  region_names <- c("Core","NBP","Functional_loop","Beta_sheets","Surface","Alpha_helices")
  for (region_name in region_names) {region_residues <- switch(
      region_name,Core = c( 4,6,7,8,9,10,11,14,15,16,17,18,19,20,21,22,23,24,40,
        42,44,46,51,52,53,54,55,56,57,58,68,71,72,75,77,78,79,
        80,81,82,83,84,89,90,92,93,96,97,99,100,101,103,109,110,
        111,112,113,114,115,116,118,125,130,133,134,137,139,141,
        142,143,144,145,146,151,152,155,156,157,158,159,160,162,163),
      NBP = c(12,13,14,15,16,17,18,28,29,30,32,34,35,57,60,61,116,117,119,120,145,146,147),
      Functional_loop = c(10,11,12,13,14,15,16,17,## p-loop
        25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,## switch i
        58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76),## switch ii
      Beta_sheets = c(3:9, 38:44, 51:57, 77:84, 109:115, 139:143),
      Surface = c( 1,2,3,5,12,13,25:39,41,43,45,47:50,59:67,
        69,70,73,74,76,85:88,91,94,95,98,102,104:108,
        117,119:124,126:129,131,132,135,136,138,140,
        147:150,153,154,161,164:188),
      Alpha_helices = c(15:24, 67:73, 87:104, 127:136, 148:166))
    p <- multimodalallostery_plot_protein_pair_region(analysis_result = result,region_residues = region_residues,region_name = region_name,
      point_size = 2.5,alpha = 0.7,base_size = 12,xlim = c(-1.5, 3),ylim = c(-1.5, 3))
    p_with_legend <- p +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(color = ggplot2::guide_legend(ncol = 2,byrow = TRUE))
    save_pdf( p_with_legend,"F6b",region_name,4.1,5.5)}}

run_F6b()



# Supplementary Figure 1: library agreement, replicate agreement, heatmaps.
## SF1b
run_SF1b <- function() {
  dir.create(file.path(OUTPUT_DIR, "SupplementaryFigure1"), recursive = TRUE,
             showWarnings = FALSE)
  comparisons <- list(
    Abundance = c("CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
      "CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
      "MA_RAS_abundance_2_fitness_replicates_fullseq.RData",
      "MA_RAS_abundance_3_fitness_replicates_fullseq.RData"),
    RAF1 = c("CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
      "CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData",
      "MA_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
      "MA_RAS_binding_RAF_3_fitness_replicates_fullseq.RData"),
    K55 = c("CW_RAS_binding_K55_2_fitness_replicates_fullseq.RData",
      "CW_RAS_binding_K55_3_fitness_replicates_fullseq.RData",
      "MA_RAS_binding_K55_2_fitness_replicates_fullseq.RData",
      "MA_RAS_binding_K55_3_fitness_replicates_fullseq.RData"),
    K27 = c("CW_RAS_binding_K27_2_fitness_replicates_fullseq.RData",
      "CW_RAS_binding_K27_3_fitness_replicates_fullseq.RData",
      "MA_RAS_binding_K27_2_fitness_replicates_fullseq.RData",
      "MA_RAS_binding_K27_3_fitness_replicates_fullseq.RData"))
  for (assay in names(comparisons)) {
    files <- file.path(FITNESS_DIR, comparisons[[assay]])
    invisible(lapply(files, require_file, label = paste0("SF1b ", assay, " fitness")))
    multimodalallostery_compare_fitness_libraries_singlemut_overall_no_block1(
      files[1], files[2], files[3], files[4], wt_aa,
      output_file = file.path(OUTPUT_DIR, "SupplementaryFigure1", paste0("SF1b_", assay, ".pdf")),
      x_lab = paste(assay, "nicking library fitness"),
      y_lab = paste(assay, "synthetic library fitness"),
      main_title = "Comparison of fitness data between synthetic and nicking libraries\nsingle mutation")
  }
}

run_SF1b()

## SF1c
sf1_fitness_blocks <- list(
  K13 = fitness_blocks$K13, K19 = fitness_blocks$K19,
  stability = fitness_blocks$Abundance, RAF1 = fitness_blocks$RAF1,
  K55 = file.path(FITNESS_DIR, c("CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData",
    "MA_RAS_binding_K55_2_fitness_replicates_fullseq.RData", "MA_RAS_binding_K55_3_fitness_replicates_fullseq.RData")),
  K27 = file.path(FITNESS_DIR, c("CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData",
    "MA_RAS_binding_K27_2_fitness_replicates_fullseq.RData", "MA_RAS_binding_K27_3_fitness_replicates_fullseq.RData")))

run_SF1c <- function() {
  colour_scheme <- list(blue = "#1B38A6", red = "#F4270C", orange = "#F4AD0C",
    green = "#09B636", yellow = "#F1DD10", purple = "#C68EFD",
    `hot pink` = "#FF0066", `light blue` = "#75C2F6", `light red` = "#FF6A56",
    `dark red` = "#A31300", `dark green` = "#007A20", pink = "#FFB0A5")
  for (assay in names(sf1_fitness_blocks)) {
    files <- sf1_fitness_blocks[[assay]]
    invisible(lapply(files, require_file, label = paste0("SF1c ", assay, " fitness")))
    normalized <- krasddpcams::krasddpcams__normalize_growthrate_fitness(
      block1_dimsum_df = files[1], block2_dimsum_df = files[2], block3_dimsum_df = files[3])
    p <- if (assay %in% c("K13", "K19"))
      multimodalallostery_plot_fitness_correlation_blocks_BI2(normalized, assay, colour_scheme)
    else multimodalallostery_plot_fitness_correlation_blocks_BI1(normalized, assay, colour_scheme)
    save_pdf(p, "SF1c", assay, 4, 4)
  }
}

run_SF1c()

## SF1d
run_SF1d <- function() {
  titles <- c(stability = "KRAS-Abundance", RAF1 = "KRAS-RAF1", K13 = "KRAS-DARPin K13",
              K19 = "KRAS-DARPin K19", K55 = "KRAS-DARPin K55", K27 = "KRAS-DARPin K27")
  for (assay in names(titles)) {
    files <- sf1_fitness_blocks[[assay]]
    invisible(lapply(files, require_file, label = paste0("SF1d ", assay, " fitness")))
    normalized <- wlab.block::nor_fitness(block1 = files[1], block2 = files[2], block3 = files[3])
    single <- wlab.block::nor_fitness_single_mut(normalized)
    single <- wlab.block::pos_id(single, wt_aa)
    save_pdf(multimodalallostery_fitness_heatmap(single, wt_aa, titles[[assay]], c(-2, 1.5)),
             "SF1d", assay, 20, 6)
  }
}

run_SF1d()

# Supplementary Figure 2: MoCHI model evaluation and external-energy checks.
sf2_abundance_blocks <- file.path(FITNESS_DIR, c(
  "CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
  "CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
  "MA_RAS_abundance_2_fitness_replicates_fullseq.RData",
  "CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
  "MA_RAS_abundance_3_fitness_replicates_fullseq.RData"))

sf2_folding_data <- function() {
  invisible(lapply(sf2_abundance_blocks, require_file, label = "SF2 abundance fitness"))
  multimodalallostery_merge_ddgf_fitness_blocks(
    prediction = require_file(PREDICTIONS_FILE, "MoCHI predictions"),
    folding_ddG = folding_weight(), block1_dimsum_df = sf2_abundance_blocks[1],
    block2_dimsum_df = sf2_abundance_blocks[2], block3_dimsum_df = sf2_abundance_blocks[3],
    block4_dimsum_df = sf2_abundance_blocks[4], block5_dimsum_df = sf2_abundance_blocks[5],
    wt_aa_input = wt_aa)
}

## SF2a
run_SF2a <- function() {
  dat <- sf2_folding_data()
  linear_files <- file.path(WEIGHTS_DIR, c("linears_weights_Abundance1.txt",
    "linears_weights_Abundance2_1.txt", "linears_weights_Abundance2_2.txt",
    "linears_weights_Abundance3_1.txt", "linears_weights_Abundance3_2.txt"))
  for (i in seq_along(linear_files)) {
    require_file(linear_files[i], "SF2a linear weights")
    p <- krasddpcams::krasddpcams__plot2d_ddGf_fitness(
      pre_nor = dat, fold_n = 1, mochi_parameters = linear_files[i], phenotypen = i,
      RT = 0.001987 * (273 + 30), bin_input = 50)
    save_pdf(p, "SF2a", paste0("block", i), 60 / 25.4, 35 / 25.4)
  }
}


run_SF2a()


## SF2b
run_SF2b <- function() {
  dat <- sf2_folding_data()
  for (i in 1:5) save_pdf(multimodalallostery_plot_ddg_fitness_per_block(dat, i, TRUE),
                           "SF2b", paste0("block", i), 60 / 25.4, 45 / 25.4)
}

run_SF2b()

## SF2c/SF2d
run_sf2_binding <- function(assay, panel, phenotypes) {
  files <- fitness_blocks[[assay]]
  invisible(lapply(files, require_file, label = paste0(panel, " ", assay, " fitness")))
  dat <- multimodalallostery_k13_19_get_ob_pre_fitness_binding_correlation_3blocks(
    prediction = require_file(PREDICTIONS_FILE, "MoCHI predictions"),
    block1_dimsum_df = files[1], block2_dimsum_df = files[2], block3_dimsum_df = files[3],
    assay_sele = assay, wt_aa_input = wt_aa)
  for (i in seq_along(phenotypes))
    save_pdf(multimodalallostery_plot_ddg_fitness_per_block(dat, phenotypes[i], TRUE),
             panel, paste0("block", i), 60 / 25.4, 45 / 25.4)
}

run_SF2c <- function() run_sf2_binding("K13", "SF2c", 6:8)
run_SF2c()

run_SF2d <- function() run_sf2_binding("K19", "SF2d", 9:11)
run_SF2d()

## SF2e
run_SF2e <- function() {
  p <- krasddpcams::krasddpcams__plot_RAF_invitro_cor(
    input = weight("RAF1"), assay_name = "RAF1",
    colour_scheme <- list(red = "#F4270C",blue = "#1B38A6"))
  save_pdf(p, "SF2e", "RAF1_invitro", 50 / 25.4, 50 / 25.4)
}

run_SF2e()

## SF2f
run_SF2f <- function() {
  assays <- c("K55", "K27")
  for (assay in assays) {
    path <- if (assay == "Folding") folding_weight() else weight(assay)
    save_pdf(multimodalallostery_dd_g_heatmap( path, wt_aa, title = assay), "SF2f", assay,
             width = 20, height = 6)
  }
}

run_SF2f()

## SF2g

run_SF2g <- function() {
  current <- c(Folding = folding_weight(), RAF1 = weight("RAF1"),K55 = weight("K55"),K27 = weight("K27"))
  reference <- c(Folding = file.path(WEIGHTS_DIR, "weights_Folding_weng.txt"),
                 RAF1 = file.path(WEIGHTS_DIR, "weights_Binding_RAF_weng.txt"),
                 K55 = file.path(WEIGHTS_DIR, "weights_Binding_K55_weng.txt"),
                 K27 = file.path(WEIGHTS_DIR, "weights_Binding_K27_weng.txt"))
  for (assay in names(current)) {current_file <- current[[assay]]
                                 reference_file <- reference[[assay]]
    require_file(current_file, paste0("SF2g current ", assay, " weights"))
    require_file(reference_file, paste0("SF2g Weng ", assay, " weights"))
    x <- data.table::fread(current_file)[, c(1, 3, 20:22)]
    y <- data.table::fread(reference_file)[, c(1, 3, 20:22)]
    p <- multimodalallostery_plot_dd_g_correlation(x, y, paste0(assay, "_this study"), paste0(assay, "_Weng"),limits = c(-1.6, 2.8))
    save_pdf(p, "SF2g", assay, 4, 4)}
}

run_SF2g()

# Supplementary Figure 3: enrichment and alternate assay/spatial views.
## SF3a
run_SF3a <- function() {
  for (assay in c("RAF1", "K13", "K19")) {
    result <- multimodalallostery_process_single_assay( weight(assay), assay,
                       input("annotation_5"), input("contact_shell"), threshold = 0.40)
    for (part in c("NBP_results", "second_shell_results", "beta_sheet_results"))
      if (!is.null(result[[part]])) save_csv(result[[part]], "SF3a", paste0(assay, "_", part))
    prepared <- multimodalallostery_prepare_plot_data( result$NBP_results,
                         result$second_shell_results, result$beta_sheet_results, assay)
    out <- file.path(OUTPUT_DIR, "SupplementaryFigure3")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    multimodalallostery_create_enrichment_plot( prepared, assay, out)
  }
}

run_SF3a()

## SF3b

run_SF3b <- function() {
  dat <- multimodalallostery_load_dd_g_data(weight("RAF1"), "RAF1")
  beta <- multimodalallostery_prepare_beta_sheet_ddg(
    dat, data.table::fread(input("annotation_5")), sheet_rects)
  save_pdf(multimodalallostery_plot_beta_sheet_ddg(beta, "Binding free energy change(RAF1) \n(kcal/mol)"),
    "SF3b", "RAF1", 4, 4)
}

run_SF3b()

## SF3c
run_SF3c <- function() {
  result <- multimodalallostery_manhatta_plot_single_assay( weight("K19"), "K19",
                     input("annotation_5"), sheet_rects, helix_rects, wt_aa,
                     threshold = 0.4)
  save_pdf(result$plot, "SF3c", "K19", 16, 6)
}

run_SF3c()


## SF3e
run_SF3e <- function() save_pdf(multimodalallostery_plot_energy_distance_decay_expfit(
                                     weight("K19"), "K19", input("annotation_5")), "SF3e", "K19")
run_SF3e()

## SF3f
run_SF3f <- function() save_pdf(multimodalallostery_plot_energy_distance_decay_expfit_directional(
                                     weight("K19"), "K19", input("annotation_5")), "SF3f", "K19")
run_SF3f()

## SF3g
run_SF3g <- function() save_pdf(multimodalallostery_plot_energy_distance_decay_expfit_contact_shell(
                                     weight("K19"), "K19", data.table::fread(input("contact_shell"))), "SF3g", "K19")
run_SF3g()

## SF3h
run_SF3h <- function() save_pdf(multimodalallostery_plot_energy_distance_decay_expfit_contact_shell_directional(
                                     weight("K19"), "K19", data.table::fread(input("contact_shell"))), "SF3h", "K19")
run_SF3h()

# Supplementary Figure 4: weighted mean over eight assays, with display subset.
## SF4a
run_SF4a <- function() {
  assays <- c("RAF1", "K13", "K19", "K27", "K55", "RAL", "PI3", "SOS")
  files <- stats::setNames(vapply(assays, weight, character(1)), assays)
  result <- multimodalallostery_plot_weighted_mean_dd_g_distance_with_cross( files,
                     all_assays = assays, plot_assays = c("RAF1", "K13", "K19"),
                     anno_file = input("annotation_8"))
  save_pdf(result$plot, "SF4a", "weighted_distance", 10, 6)
  save_csv(result$data, "SF4a", "weighted_distance_data")
}

run_SF4a()

# Supplementary Figure 5: pairwise correlations and cross-assay mapping.
## SF5a
run_SF5a <- function() {
  assays <- c("RAF1", "RAL", "PI3", "SOS","K55", "K27", "K13", "K19")
  files <- stats::setNames(lapply(assays, weight), assays)
  result <- multimodalallostery_calculate_correlation_matrix( files, assays)
  heatmap <- multimodalallostery_plot_correlation_heatmap( result$cor_matrix, result$p_matrix)
  save_pdf(heatmap$gtable,
           "SF5a", "correlation_matrix", 10, 8)
}
run_SF5a()

## SF5b
run_SF5b <- function() {
  path <- require_file(file.path(MANUSCRIPT_FIGURE_DIR,"mutation_classification_summary_8binder_fixed_threshold_0.40.csv"), "SF5b pair counts")
  pair_order <- c(
    "RAF1 vs K13", "RAF1 vs K19",
    "RALGDS vs K13", "RALGDS vs K19",
    "PI3KCG vs K13", "PI3KCG vs K19",
    "SOS1 vs K13", "SOS1 vs K19",
    "K55 vs K13", "K55 vs K19",
    "K27 vs K13", "K27 vs K19")
  prepared <- multimodalallostery_prepare_allosteric_pair_counts(data.table::fread(path),pair_order)
  save_pdf(multimodalallostery_plot_allosteric_pair_counts(prepared,c("#1B38A6", "#75C2F6")),"SF5b","pair_counts",5,3)
}

run_SF5b()

## SF5c
run_SF5c <- function() {
  anti <- multimodalallostery_get_raf1_vs_k27_anticorrelated( weight("RAF"), weight("K27"),
                   input("annotation_8"), fixed_threshold = 0.4)
  mapped <- multimodalallostery_load_k13_vs_k19_data( weight("K13"), weight("K19"),
                     input("annotation_8"))
  save_pdf(multimodalallostery_plot_mapped_anticorrelated( mapped, anti), "SF5c", "mapped")
  save_csv(anti, "SF5c", "anticorrelated_mutations")
}

run_SF5c()

## SF5d
run_SF5d_left <- function() {
  base <- file.path(MANUSCRIPT_FIGURE_DIR)
  paths <- file.path(base, c(
    "RAF1_vs_K27_anticorrelated_mutations_fixed_threshold_0.40.csv",
    "RAF1_vs_K13_anticorrelated_mutations_fixed_threshold_0.40.csv",
    "K27_vs_K13_anticorrelated_mutations_fixed_threshold_0.40.csv"))
  invisible(lapply(paths, require_file, label = "SF5d pairwise mutation table"))
  tables <- stats::setNames(lapply(paths, data.table::fread),
    c("RAF1 vs K27", "RAF1 vs K13", "K27 vs K13"))
  sets <- multimodalallostery_prepare_mutation_sets(tables)
  save_pdf(multimodalallostery_plot_mutation_venn(sets), "SF5d", "venn_left")
}
run_SF5d_left()

run_SF5d_right <- function() {
  base <- file.path(MANUSCRIPT_FIGURE_DIR)
  paths <- file.path(base, c(
    "RAF1_vs_K27_anticorrelated_mutations_fixed_threshold_0.40.csv",
    "RAF1_vs_K19_anticorrelated_mutations_fixed_threshold_0.40.csv",
    "K27_vs_K19_anticorrelated_mutations_fixed_threshold_0.40.csv"))
  invisible(lapply(paths, require_file, label = "SF5d pairwise mutation table"))
  tables <- stats::setNames(lapply(paths, data.table::fread),
                            c("RAF1 vs K27", "RAF1 vs K19", "K27 vs K19"))
  sets <- multimodalallostery_prepare_mutation_sets(tables)
  save_pdf(multimodalallostery_plot_mutation_venn(sets), "SF5d", "venn_right")
}
run_SF5d_right()


## SF5e
run_SF5e <- function() {
  both <- multimodalallostery_get_k13_vs_k19_both_promoting( weight("K13"), weight("K19"),
                   input("annotation_8"), fixed_threshold = 0.4)
  x <- multimodalallostery_load_k27_vs_k13_data( weight("K27"), weight("K13"), input("annotation_8"))
  y <- multimodalallostery_load_k27_vs_k19_data( weight("K27"), weight("K19"), input("annotation_8"))
  save_pdf(multimodalallostery_plot_mapped_both_promoting_k27_k13( x, both), "SF5e", "K27_K13")
  save_pdf(multimodalallostery_plot_mapped_both_promoting_k27_k19( y, both), "SF5e", "K27_K19")
  save_csv(both, "SF5e", "both_promoting_mutations")
}

run_SF5e()

# Supplementary Figure 6: region enrichment, residue correlations and site plots.
## SF6a
run_SF6a <- function() {
  assays <- c("RAF1", "K55", "K27", "K13", "K19")
  files <- stats::setNames(lapply(assays, function(x) weight(if (x == "RAF1") "RAF" else x)), assays)
  regions <- defs$region_sets
  result <- multimodalallostery_run_pair_region_enrichment( "RAF1 vs K13", files,
                     input("annotation_8"), regions, fixed_threshold = 0.4)
  prepared <- multimodalallostery_prepare_region_enrichment_plot_data( result, "RAF1 vs K13",
                       names(regions))
  save_pdf(multimodalallostery_plot_region_enrichment( prepared,
                    c(Correlated = "#F4AD0C", `Anti-correlated` = "#1B38A6",
                      Other = "grey80")), "SF6a", "region_enrichment", 10, 6)
  save_csv(result$or, "SF6a", "odds_ratios")
  save_csv(result$plot, "SF6a", "fractions")
}

run_SF6a()

## SF6b
run_SF6b <- function() {
  assays <- c("RAF1", "K13")
  files <- stats::setNames(lapply(assays, function(x) weight(if (x == "RAF1") "RAF" else x)), assays)
  regions <- defs$region_sets
  binding_sites_map <- list(
    RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
    K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99,
            101, 102, 105, 106, 107, 129, 133, 136, 137, 138))
  contact_shell <- data.table::fread(input("contact_shell"))
  second_shell_map <- multimodalallostery_prepare_second_shell_map(contact_shell, assays)
  result <- multimodalallostery_run_pair_second_shell_enrichment("RAF1 vs K13", files,input("annotation_8"), regions,
    second_shell_map, binding_sites_map,fixed_threshold = 0.4)
  save_pdf(result$plot, "SF6b", "second_shell_enrichment", 10, 6)
  save_csv(result$or, "SF6b", "odds_ratios")
  save_csv(result$fractions, "SF6b", "fractions")
}

run_SF6b()


## SF6c
run_SF6c <- function() {
  path <- require_file(file.path(MANUSCRIPT_FIGURE_DIR,
    "RAF1_vs_K13_correlation_results_with_WT_0.csv"))
  data <- data.table::fread(path)
  prepared <- multimodalallostery_prepare_residue_correlations(data)
  tracks <- multimodalallostery_prepare_annotation_tracks(defs$annotation_tracks, prepared)
  p <- multimodalallostery_combine_residue_panels(
    multimodalallostery_plot_residue_correlations(prepared),
    multimodalallostery_plot_annotation_tracks(tracks, levels(prepared$residue_order)))
  save_pdf(p, "SF6c", "residue_correlations", 16, 6)
}

run_SF6c()

## SF6d
run_SF6d <- function() {
  result <- multimodalallostery_analyze_all_sites( weight("RAF"), weight("K13"),
                     "RAF1", "K13", input("annotation_8"), fixed_threshold = 0.4)
  for (site in c(145, 15, 48, 55, 77, 163))
    save_pdf(multimodalallostery_plot_site_scatter( result, target_position = site,
                      point_size = 5, base_size = 20, xlim = c(-1.5, 3),
                      ylim = c(-1.5, 3)), "SF6d", paste0("site_", site))
}

run_SF6d()

# ---- Run the current 33 panels in the required order -------------------------
run_all_figures <- function(panels = NULL) {
  runners <- list(
    F1e = run_F1e, F1f = run_F1f, F1g = run_F1g, F1h = run_F1h, F1i = run_F1i,
    F2a = run_F2a, F2c = run_F2c, F2d = run_F2d,
    F3a = run_F3a, F3c = run_F3c, F3d = run_F3d, F3e = run_F3e, F3f = run_F3f,
    F4b = run_F4b, F5a = run_F5a, F6b = run_F6b,
    SF1b = run_SF1b, SF1c = run_SF1c, SF1d = run_SF1d,
    SF2a = run_SF2a, SF2b = run_SF2b, SF2c = run_SF2c, SF2d = run_SF2d,
    SF2e = run_SF2e, SF2g = run_SF2g,
    SF3a = run_SF3a, SF3b = run_SF3b, SF3c = run_SF3c, SF3e = run_SF3e,
    SF3f = run_SF3f, SF3g = run_SF3g, SF3h = run_SF3h,
    SF4a = run_SF4a,
    SF5a = run_SF5a, SF5b = run_SF5b, SF5c = run_SF5c, SF5d = run_SF5d, SF5e = run_SF5e,
    SF6a = run_SF6a, SF6b = run_SF6b, SF6c = run_SF6c, SF6d = run_SF6d)
  if (!is.null(panels)) {
    unknown <- setdiff(panels, names(runners))
    if (length(unknown)) stop("Unknown panel: ", paste(unknown, collapse = ", "))
    runners <- runners[panels]
  }
  if (!length(runners)) stop("No panels selected.")
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  for (panel in names(runners)) {
    message("[", panel, "] starting")
    runners[[panel]]()
    message("[", panel, "] completed")
  }
  message("Workflow completed: ", length(runners), " panel(s). Output: ", OUTPUT_DIR)
  invisible(names(runners))
}

run_all_figures(PANELS_TO_RUN)
