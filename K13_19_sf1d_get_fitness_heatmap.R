# Fitness Heatmap Visualization for KRAS Mutations
# This script creates heatmaps to visualize fitness effects of single amino acid mutations

library(wlab.block)
library(data.table)
library(ggplot2)

#' Create fitness heatmap for single amino acid mutations
#'
#' @param input Normalized fitness data for single mutations
#' @param wt_aa Wild-type amino acid sequence
#' @param title Plot title
#' @param legend_limits Limits for color legend
#' @return ggplot object showing fitness heatmap

nor_fitness_5blocks <- function(required_file) {
  if(length(required_file) != 5) stop("必须提供 5 个文件路径")
  
  all_data <- list()
  gr_lm <- data.frame(matrix(nrow = 2, ncol = 0))
  fit_lm <- data.frame(matrix(nrow = 2, ncol = 0))
  
  # 遍历 5 个 block
  for (i in 1:5) {
    load(required_file[i])
    
    # 加权平均 growthrate
    all_variants$gr_over_sigmasquared <- all_variants$growthrate / (all_variants$growthrate_sigma)^2
    all_variants$one_over_sigmasquared <- 1 / (all_variants$growthrate_sigma)^2
    stop1 <- sum(all_variants[STOP == TRUE, ]$gr_over_sigmasquared, na.rm = TRUE) / 
      sum(all_variants[STOP == TRUE, ]$one_over_sigmasquared, na.rm = TRUE)
    wt1 <- sum(all_variants[WT == TRUE, ]$gr_over_sigmasquared, na.rm = TRUE) / 
      sum(all_variants[WT == TRUE, ]$one_over_sigmasquared, na.rm = TRUE)
    gr_lm <- dplyr::bind_cols(gr_lm, c(stop = stop1, wt = wt1))
    colnames(gr_lm)[i] <- paste0("block", i)
    
    # 加权平均 fitness
    all_variants$fitness_over_sigmasquared <- all_variants$fitness / (all_variants$sigma)^2
    all_variants$one_over_fitness_sigmasquared <- 1 / (all_variants$sigma)^2
    stop1_fitness <- sum(all_variants[STOP == TRUE, ]$fitness_over_sigmasquared, na.rm = TRUE) / 
      sum(all_variants[STOP == TRUE, ]$one_over_fitness_sigmasquared, na.rm = TRUE)
    wt1_fitness <- sum(all_variants[WT == TRUE, ]$fitness_over_sigmasquared, na.rm = TRUE) / 
      sum(all_variants[WT == TRUE, ]$one_over_fitness_sigmasquared, na.rm = TRUE)
    fit_lm <- dplyr::bind_cols(fit_lm, c(stop = stop1_fitness, wt = wt1_fitness))
    colnames(fit_lm)[i] <- paste0("block", i)
    
    all_data[[i]] <- all_variants
  }
  
  # 给第一个 block 添加 nor 列
  all_data[[1]][, `:=`(nor_gr, growthrate)]
  all_data[[1]][, `:=`(nor_gr_sigma, growthrate_sigma)]
  all_data[[1]][, `:=`(nor_fitness, fitness)]
  all_data[[1]][, `:=`(nor_fitness_sigma, sigma)]
  
  # 对后续 block 做归一化
  for (i in 2:5) {
    # growthrate
    formula1 <- as.formula(paste0("block1 ~ block", i))
    b1b2 <- summary(lm(formula = formula1, data = gr_lm))
    a2 <- b1b2$coefficients[2,1]
    b2 <- b1b2$coefficients[1,1]
    
    # fitness
    formula2 <- as.formula(paste0("block1 ~ block", i))
    c1c2 <- summary(lm(formula = formula2, data = fit_lm))
    d2 <- c1c2$coefficients[2,1]
    e2 <- c1c2$coefficients[1,1]
    
    all_data[[i]][, `:=`(nor_gr, growthrate * a2 + b2)]
    all_data[[i]][, `:=`(nor_gr_sigma, growthrate_sigma * a2)]
    all_data[[i]][, `:=`(nor_fitness, fitness * d2 + e2)]
    all_data[[i]][, `:=`(nor_fitness_sigma, sigma * d2)]
  }
  
  # 合并所有 block
  data_after_nor <- dplyr::bind_rows(all_data, .id = "block")
  return(data_after_nor)
}



fitness_heatmap <- function(input, wt_aa, title = "fitness", legend_limits = c(-1.5, 1)) {
  # 定义氨基酸顺序
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  num <- nchar(wt_aa) + 1
  
  # 准备单点突变数据
  input_single <- input
  input_single[, position := AA_Pos1]
  input_single[, WT_AA := wtcodon1]
  
  # 构建热图模板
  heatmap_tool_fitness <- data.table(
    wtcodon1 = rep(unlist(strsplit(wt_aa, "")), each = 21),
    position = rep(2:num, each = 21),
    codon1 = c(unlist(aa_list), "*")
  )
  
  # 合并实验数据
  heatmap_tool_fitness_anno_single <- merge(
    input_single, heatmap_tool_fitness,
    by = c("wtcodon1", "position", "codon1"), all = TRUE
  )
  
  # 设置氨基酸顺序
  heatmap_tool_fitness_anno_single <- within(
    heatmap_tool_fitness_anno_single,
    codon1 <- factor(codon1, levels = c(
      "*", "D", "E", "R", "H", "K", "S", "T", "N", "Q",
      "C", "G", "P", "A", "V", "I", "L", "M", "F", "W", "Y"
    ))
  )
  
  # WT fitness = 0
  heatmap_tool_fitness_anno_single[wtcodon1 == codon1, nor_fitness_nooverlap := 0]
  
  # 绘图
  ggplot() +
    theme_classic() +
    geom_tile(
      data = heatmap_tool_fitness_anno_single[position > 1, ],
      aes(x = position, y = codon1, fill = nor_fitness_nooverlap)
    ) +
    scale_x_discrete(limits = c(2:num), labels = c(2:num)) +
    theme(
      axis.text.x = element_text(
        size = 8, vjust = 0.5, hjust = 0.5,
        color = c(NA, NA, NA, rep(c("black", NA, NA, NA, NA), nchar(wt_aa) %/% 5))
      )
    ) +
    scale_fill_gradient2(
      limits = legend_limits,
      low = "#F4270C", mid = "gray", high = "#1B38A6",
      name = "Fitness",
      midpoint = 0,
      na.value = "white",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    ylab("Mutant aa") +
    ggtitle(title) +
    labs(fill = NULL) +
    geom_text(
      data = heatmap_tool_fitness_anno_single[position > 1 & wtcodon1 == codon1, ],
      aes(x = position, y = codon1),
      label = "-", size = 3
    ) +
    theme(
      text = element_text(size = 8),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = c(1, 1.38),
      title = element_text(size = 8),
      legend.justification = c(1, 1),
      legend.direction = "horizontal",
      legend.text = element_text(size = 8),
      axis.title.x = element_text(size = 8, face = "plain"),
      axis.title.y = element_text(size = 8, face = "plain"),
      axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(
        family = "Courier", angle = 90, size = 9.5,
        vjust = 0.5, hjust = 0.5,
        margin = margin(0, -0.5, 0, 0, "mm")
      ),
      legend.key.height = unit(3.1, "mm"),
      legend.key.width = unit(4, "mm"),
      legend.key.size = unit(1, "mm"),
      plot.margin = margin(0, -0, 0, 0)
    ) +
    coord_fixed()
}

#--------------------
#调用函数
#--------------------

wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
####Abundance
nor_fit<-nor_fitness(block1="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                     block2="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/Abundance_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                     block3="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/Abundance_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData")

nor_fit_single<-nor_fitness_single_mut(input=nor_fit)
nor_fit_single<-pos_id(nor_fit_single,wt_aa)
#step 3
fitness_heatmap(nor_fit_single,wt_aa,title = "KRAS-Abundance" ,legend_limits = c(-2, 1.5))


ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figures1/KRAS-Abundance.pdf",height = 6,width=20)








####RAF1
nor_fit<-nor_fitness(block1="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
                     block2="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/RAF_block2_Q20_rbg_filter2_20250829_fitness_replicates.RData",
                     block3="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/RAF_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")

nor_fit_single<-nor_fitness_single_mut(input=nor_fit)
nor_fit_single<-pos_id(nor_fit_single,wt_aa)
#step 3
fitness_heatmap(nor_fit_single,wt_aa,title = "KRAS-RAF1",legend_limits = c(-2, 1.5))


ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figures1/KRAS-RAF1.pdf",height = 6,width=20)









####K13
nor_fit<-nor_fitness(block1="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K13_block1_Q20_rbg_filter2_20251109_fitness_replicates.RData",
                     block2="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K13_block2_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData",
                     block3="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K13_block3_Q20_rbg_filter2_20250829_fitness_replicates_cleaned.RData")

nor_fit_single<-nor_fitness_single_mut(input=nor_fit)
nor_fit_single<-pos_id(nor_fit_single,wt_aa)
#step 3
fitness_heatmap(nor_fit_single,wt_aa,title = "KRAS-DARPin K13",legend_limits = c(-2, 1.5))


ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figures1/KRAS-DARPin K13.pdf",height = 6,width=20)





####K19
nor_fit<-nor_fitness(block1="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K19_block1_Q20_rbg_filter8_20251109_fitness_replicates_cleaned.RData",
                     block2="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K19_block2_Q20_rbg_filter1_20251107_fitness_replicates_cleaned.RData",
                     block3="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K19_block3_Q20_rbg_filter2_20250829_fitness_replicates.RData")

nor_fit_single<-nor_fitness_single_mut(input=nor_fit)
nor_fit_single<-pos_id(nor_fit_single,wt_aa)
#step 3
fitness_heatmap(nor_fit_single,wt_aa,title = "KRAS-DARPin K19",legend_limits = c(-2, 1.5))


ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figures1/KRAS-DARPin K19.pdf",height = 6,width=20)





####K55
nor_fit<-nor_fitness(block1="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData",
                     block2="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K55_block2_Q20_rbg2_filter1_fitness_replicates.RData",
                     block3="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K55_block3_Q20_rbg_1_filter1_fitness_replicates.RData")

nor_fit_single<-nor_fitness_single_mut(input=nor_fit)
nor_fit_single<-pos_id(nor_fit_single,wt_aa)
#step 3
fitness_heatmap(nor_fit_single,wt_aa,title = "KRAS-DARPin K55",legend_limits = c(-2, 1.5))


ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figures1/KRAS-DARPin K55.pdf",height = 6,width=20)






####K27
nor_fit<-nor_fitness(block1="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData",
                     block2="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K27_block2_Q20_rbg3_filter1_fitness_replicates.RData",
                     block3="C:/Users/36146/OneDrive - USTC/DryLab/fitness RData/K27_block3_Q20_rbg3_filter1_fitness_replicates.RData")

nor_fit_single<-nor_fitness_single_mut(input=nor_fit)
nor_fit_single<-pos_id(nor_fit_single,wt_aa)
#step 3
fitness_heatmap(nor_fit_single,wt_aa,title = "KRAS-DARPin K27",legend_limits = c(-2, 1.5))


ggplot2::ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260110_version/figures1/KRAS-DARPin K27.pdf",height = 6,width=20)
