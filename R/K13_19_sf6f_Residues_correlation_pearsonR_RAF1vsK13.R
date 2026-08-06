# ============================================================================
# Step 0. Packages
# ============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


# ============================================================================
# Step 1. Read data
# ============================================================================

R_value_pairwise <- fread(
  "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result5_figure6/RAF1_vs_K13_correlation_results_with_WT_0.csv"
)


# ============================================================================
# Step 2. Structural annotations
# ============================================================================

annotation_tracks <- list(
  
  BI1 = c(
    5, 21, 24, 25, 27, 29, 31, 33, 36, 37, 38, 39, 40, 41, 43,
    52, 54, 56, 64, 66, 67, 70, 71, 73, 74
  ),
  
  BI2 = c(
    63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99,
    101, 102, 105, 106, 107, 108,
    125, 129, 133, 136, 137, 138
  ),
  `β-strand 1` = 3:9,
  `β-strand 2` = 38:44,
  `β-strand 3` = 51:57,
  `β-strand 4` = 77:84,
  `β-strand 5` = 109:115,
  `β-strand 6` = 139:143,
  NBP = c(
    12, 13, 14, 15, 16, 17, 18,28, 29, 30, 32, 34, 35,57, 60, 61,116, 117, 119, 120,145, 146, 147),
  `γ-phosphate contact` = c(16, 35, 60, 12, 13, 34))
# ============================================================================
# Step 3. Prepare Pearson R data
# ============================================================================
plot_df <- R_value_pairwise %>%
  arrange(R) %>%
  mutate(
    residue_order = factor(Pos_real,levels = Pos_real),
     sig_group = case_when(
      pvalue < 0.05 & R > 0 ~ "Positive",
      pvalue < 0.05 & R < 0 ~ "Negative",
      TRUE ~ "NS"))
# ============================================================================
# Step 4. Pearson R panel
# ============================================================================
p_R <- ggplot(plot_df,
  aes(x = residue_order,y = R)) +
  geom_hline(yintercept = 0,linetype = 2,linewidth = 0.4,colour = "grey40") +
  geom_point(aes(fill = sig_group),shape = 21,size = 2.3,
    colour = "white",stroke = 0.25) +
  scale_fill_manual(values = c(Positive = "#F4AD0C",Negative = "#75C2F6",NS = "grey80")) +
  scale_y_continuous(limits = c(-1, 1),breaks = seq(-1, 1, 0.5),expand = expansion(mult = c(0.02, 0.02))) +
  labs(y = "Residue Correlation Pearson's R between RAF1 and K13",x = NULL,fill = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    plot.margin = margin(5, 5, 0, 5))
# ============================================================================
# Step 5. Generate annotation dataframe
# ============================================================================
annotation_df <- rbindlist(
  lapply(names(annotation_tracks),
    function(track){
      data.table(
        track = track,
        Pos_real = annotation_tracks[[track]])}))
annotation_df <- merge(
  annotation_df,
  plot_df[, .(Pos_real, residue_order)],
  by = "Pos_real",
  all.x = TRUE)

annotation_df$track <- factor(
annotation_df$track,
levels = rev(names(annotation_tracks)))
annotation_df$residue_order <- factor(annotation_df$residue_order,levels = levels(plot_df$residue_order))
# ============================================================================
# Step 6. Structural annotation panel
# ============================================================================
p_annotation <- ggplot(annotation_df,
  aes(x = residue_order,y = track)) +
  geom_tile(
    width = 0.95,
    height = 0.75,
    fill = "grey70") +
  scale_x_discrete(
    drop = FALSE,
    breaks = levels(plot_df$residue_order),
    labels = levels(plot_df$residue_order)) +
  labs(x = "KRAS residue position (ordered by Pearson's R)",y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6),
    axis.ticks.x = element_line(linewidth = 0.25),
    axis.text.y = element_text(size = 9),
    plot.margin = margin(0,5,5,5))
# ============================================================================
# Step 7. Combine panels
# ============================================================================
Panel_C <- patchwork::wrap_plots(p_R, p_annotation,ncol = 1,heights = c(4,2.5))
# ============================================================================
# Step 8. Display
# ============================================================================
Panel_C
# ============================================================================
# Step 9. Save
# ============================================================================
ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result5_figure6/PanelC_PearsonR_annotation5.pdf",
 Panel_C,width = 16,height = 6,dpi = 600)
