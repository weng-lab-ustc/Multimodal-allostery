# ============================================================
# Step 0. Load packages
# ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(wlab.block)


# ============================================================
# Step 1. Function for final normalization
# ============================================================

normalize_fitness_WT_STOP <- function(block1, block2, block3){
  
  
  # ----------------------------------------------------------
  # 1. Block alignment
  # ----------------------------------------------------------
  
  nor_fit <- nor_fitness(
    block1 = block1,
    block2 = block2,
    block3 = block3
  )
  
  
  # ----------------------------------------------------------
  # 2. Calculate global WT and STOP references
  # ----------------------------------------------------------
  
  WT_mean <- nor_fit %>%
    filter(WT == TRUE) %>%
    summarise(
      WT = sum(
        nor_fitness / nor_fitness_sigma^2,
        na.rm = TRUE
      ) /
        sum(
          1 / nor_fitness_sigma^2,
          na.rm = TRUE
        )
    ) %>%
    pull(WT)
  
  
  STOP_mean <- nor_fit %>%
    filter(STOP == TRUE) %>%
    summarise(
      STOP = sum(
        nor_fitness / nor_fitness_sigma^2,
        na.rm = TRUE
      ) /
        sum(
          1 / nor_fitness_sigma^2,
          na.rm = TRUE
        )
    ) %>%
    pull(STOP)
  
  
  cat("WT reference =", WT_mean, "\n")
  cat("STOP reference =", STOP_mean, "\n")
  
  
  # ----------------------------------------------------------
  # 3. Normalize to WT=0 STOP=-1
  # ----------------------------------------------------------
  
  
  nor_fit <- nor_fit %>%
    mutate(
      fitness_normalized =
        (nor_fitness - WT_mean) /
        (WT_mean - STOP_mean)
    )
  
  
  return(nor_fit)
}



# ============================================================
# Step 2. Plot function (without block faceting)
# ============================================================

plot_fitness_density <- function(
    assay_type,
    block1,
    block2,
    block3,
    output_file = NULL
){
  
  
  # ----------------------------------------------------------
  # Load and normalize
  # ----------------------------------------------------------
  
  nor_fit <- normalize_fitness_WT_STOP(
    block1,
    block2,
    block3
  )
  
  
  
  # ----------------------------------------------------------
  # Mutation classification
  # ----------------------------------------------------------
  
  nor_fit_classified <- nor_fit %>%
    mutate(
      mut_type = case_when(
        
        Nham_aa == 0 &
          Nham_nt > 0 ~ "Synonymous",
        
        STOP == TRUE |
          STOP_readthrough == TRUE ~ "Stop",
        
        Nham_aa > 0 &
          indel == FALSE &
          STOP == FALSE &
          STOP_readthrough == FALSE ~ "Missense"
      )
    ) %>%
    filter(!is.na(mut_type))
  
  
  
  cat("\nMutation distribution:", assay_type,"\n")
  print(table(nor_fit_classified$mut_type))
  
  
  
  nor_fit_plot <- nor_fit_classified %>%
    mutate(
      mut_type = factor(
        mut_type,
        levels=c(
          "Synonymous",
          "Missense",
          "Stop"
        )
      )
    )
  
  
  
  # ----------------------------------------------------------
  # Overall distribution plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    nor_fit_plot,
    aes(
      x = fitness_normalized,
      color = mut_type
    )
  )+
    geom_density(
      linewidth = 1
    )+
    scale_color_manual(
      values=c(
        "Synonymous"="#09B636",
        "Missense"="#F4AD0C",
        "Stop"="#FF6A56"
      )
    )+
    labs(
      title=paste0(
        toupper(assay_type),
        " - Fitness Distribution"
      ),
      x="Normalized Fitness",
      y="Density",
      color="Mutation Type"
    )+
    xlim(-1.5, 0.5)+
    theme_classic()+
    theme(
      legend.position = "bottom",
      plot.title = element_text(
        hjust = 0.5,
        size = 12
      ),
      text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  
  
  if(!is.null(output_file)){
    
    ggsave(
      output_file,
      p,
      width = 6,
      height = 4,
      units = "in"
    )
    
    cat(
      "Saved:",
      output_file,
      "\n"
    )
  }
  
  
  return(p)
}




### Abundance

plot_fitness_density(
  block1 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
  block2 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  block3 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  assay_type = "Abundance",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result1_figure1/Abundance_normalized_density.pdf"
)



### K13
plot_fitness_density(
  block1 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block1_Q20_rbg_filter2_20251109_fitness_replicates.RData",
  block2 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  block3 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K13_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  assay_type = "K13",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result1_figure1/K13_normalized_density.pdf"
)



### K19
plot_fitness_density(
  block1 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K19_block1_Q20_rbg_filter8_20251109_fitness_replicates.RData",
  block2 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K19_block2_Q20_rbg_filter1_20251107_fitness_replicates.RData",
  block3 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/K19_block2_Q20_rbg_filter3_20250830_fitness_replicates.RData",
  assay_type = "K19",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result1_figure1/K19_normalized_density.pdf"
)




### RAF1

plot_fitness_density(
  block1 = "C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
  block2 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  block3 = "C:/Users/36146/OneDrive - USTC/DryLab/DiMSum/DiMSum_rerun_20250821/20251010_合并同义突变数据_sigma数据清洁/RAF_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData",
  assay_type = "RAF1",
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result1_figure1/RAF1_normalized_density.pdf"
)
