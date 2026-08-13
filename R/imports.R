#' Package imports
#'
#' Imports used by consolidated Figure/Panel functions.
#' @name multimodalallostery-imports
#' @import data.table
#' @import ggplot2
#' @import tidyr
#' @importFrom dplyr %>% arrange case_when count filter group_by mutate select summarise sym ungroup
#' @importFrom GGally ggpairs
#' @importFrom beeswarm beeswarm
#' @importFrom ggpubr theme_classic2
#' @importFrom ggrepel geom_text_repel
#' @importFrom pheatmap pheatmap
#' @importFrom wesanderson wes_palette
#' @importFrom wlab.block nor_fitness nor_fitness_single_mut pos_id
#' @importFrom krasddpcams krasddpcams__get_weighted_mean_abs_ddG_mutcount
#' @importFrom krasddpcams krasddpcams__nor_overlap_single_mt_fitness
#' @importFrom krasddpcams krasddpcams__pos_id krasddpcams__pvalue krasddpcams__read_ddG
#' @importFrom grDevices cairo_pdf colorRampPalette dev.off
#' @importFrom graphics abline axis boxplot legend mtext par segments text
#' @importFrom grid unit
#' @importFrom stats as.formula complete.cases cor.test fisher.test lm median nls p.adjust pnorm predict setNames sigma
#' @importFrom utils head
NULL
