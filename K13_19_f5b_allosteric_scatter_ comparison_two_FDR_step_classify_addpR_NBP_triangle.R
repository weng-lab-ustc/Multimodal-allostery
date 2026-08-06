# =========================================================
# scatter plot between pairwise comparison----General-purpose functions
# =========================================================

library(data.table)
library(ggplot2)
library(dplyr)

# =========================================================
# 全局常量定义
# =========================================================

COLOR_MAP <- c(
  "Not significant (FDR >= 0.05)" = "grey90",
  "Other (neutral in both)" = "grey90",
  "Both promoting" = "#FFB0A5",           
  "Both disrupting" = "#F4270C",          
  "Promoting in X / Disrupting in Y" = "#F4AD0C",  
  "Disrupting in X / Promoting in Y" = "#F1DD10",  
  "Allosteric only in X" = "#1B38A6",     
  "Allosteric only in Y" = "#75C2F6"      
)

LEGEND_ORDER <- c(
  "Both promoting",
  "Both disrupting",
  "Promoting in X / Disrupting in Y",
  "Disrupting in X / Promoting in Y",
  "Allosteric only in X",
  "Allosteric only in Y",
  "Other (neutral in both)",
  "Not significant (FDR >= 0.05)"
)

# 全局binding interface位点映射（可扩展）
BINDING_SITES_MAP <- list(
  RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
  K55  = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74),
  K27  = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71),
  K13  = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138),
  K19  = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
)

# NBP residues (需要标记为三角形的位点)
NBP_RESIDUES <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 
                  57, 60, 61, 116, 117, 119, 120, 145, 146, 147)

# =========================================================
# 数据读取函数
# =========================================================
load_mutation_data <- function(input, assay_sele) {
  if (is.character(input)) {
    ddG <- fread(input)
  } else {
    ddG <- as.data.table(input)
  }
  
  ddG[, Pos_real := Pos_ref + 1]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, mt := paste0(wt_codon, Pos_real, mt_codon)]
  ddG <- ddG[id != "WT"]
  
  result <- ddG[, .(mt, Pos_real, wt_codon, mt_codon, 
                    `mean_kcal/mol`, `std_kcal/mol`)]
  setnames(result, "mean_kcal/mol", "ddG")
  setnames(result, "std_kcal/mol", "ddG_std")
  result[, assay := assay_sele]
  
  return(result)
}

# =========================================================
# 计算阈值
# =========================================================
calculate_threshold <- function(data, assay_sele, anno) {
  anno_data <- fread(anno)
  
  site_ddG <- data[, .(
    mean_ddG = sum(abs(ddG) / ddG_std^2, na.rm = TRUE) / sum(1 / ddG_std^2, na.rm = TRUE)
  ), by = Pos_real]
  
  site_sigma <- data[, .(
    sigma = sqrt(1 / sum(1 / ddG_std^2, na.rm = TRUE))
  ), by = Pos_real]
  
  site_stats <- merge(site_ddG, site_sigma, by = "Pos_real")
  site_stats <- merge(site_stats, anno_data, by.x = "Pos_real", by.y = "Pos", all.x = TRUE)
  
  scHAmin_col <- paste0("scHAmin_ligand_", assay_sele)
  binding_sites <- site_stats[get(scHAmin_col) < 5, Pos_real]
  
  threshold <- site_stats[Pos_real %in% binding_sites, 
                          sum(abs(mean_ddG) / sigma^2, na.rm = TRUE) / sum(1 / sigma^2, na.rm = TRUE)]
  
  return(threshold)
}

# =========================================================
# 第一次分类：基于效应方向
# =========================================================
classify_by_direction <- function(ddG_x, ddG_y, threshold_x, threshold_y) {
  sig_x <- abs(ddG_x) > threshold_x
  sig_y <- abs(ddG_y) > threshold_y
  
  disrupt_x <- ddG_x > threshold_x
  promote_x <- ddG_x < -threshold_x
  disrupt_y <- ddG_y > threshold_y
  promote_y <- ddG_y < -threshold_y
  
  result <- rep("neutral", length(ddG_x))
  
  result[sig_x & sig_y & promote_x & promote_y] <- "Both promoting"
  result[sig_x & sig_y & disrupt_x & disrupt_y] <- "Both disrupting"
  result[sig_x & sig_y & promote_x & disrupt_y] <- "Promoting in X / Disrupting in Y"
  result[sig_x & sig_y & disrupt_x & promote_y] <- "Disrupting in X / Promoting in Y"
  result[sig_x & !sig_y] <- "Allosteric only in X"
  result[!sig_x & sig_y] <- "Allosteric only in Y"
  result[!sig_x & !sig_y] <- "Not significant (FDR >= 0.05)"
  
  return(result)
}

# =========================================================
# 第二次分类：基于FDR重新评估所有类型
# =========================================================
reclassify_by_FDR <- function(direction_class, pass_FDR_x, pass_FDR_y) {
  result <- direction_class
  
  for (i in 1:length(direction_class)) {
    
    if (direction_class[i] %in% c("Both promoting", "Both disrupting", 
                                  "Promoting in X / Disrupting in Y", 
                                  "Disrupting in X / Promoting in Y")) {
      if (!(pass_FDR_x[i] & pass_FDR_y[i])) {
        if (pass_FDR_x[i] & !pass_FDR_y[i]) {
          result[i] <- "Allosteric only in X"
        } else if (!pass_FDR_x[i] & pass_FDR_y[i]) {
          result[i] <- "Allosteric only in Y"
        } else {
          result[i] <- "Not significant (FDR >= 0.05)"
        }
      }
    }
    
    else if (direction_class[i] == "Allosteric only in X") {
      if (!pass_FDR_x[i]) {
        result[i] <- "Not significant (FDR >= 0.05)"
      }
    }
    
    else if (direction_class[i] == "Allosteric only in Y") {
      if (!pass_FDR_y[i]) {
        result[i] <- "Not significant (FDR >= 0.05)"
      }
    }
  }
  
  return(result)
}

# =========================================================
# 准备合并数据和计算双FDR
# =========================================================
prepare_merged_data_with_FDR <- function(input_x, input_y, assay_x, assay_y, anno) {
  data_x <- load_mutation_data(input_x, assay_x)
  data_y <- load_mutation_data(input_y, assay_y)
  
  threshold_x <- calculate_threshold(data_x, assay_x, anno)
  threshold_y <- calculate_threshold(data_y, assay_y, anno)
  
  # 移除binding interface位点
  if (!is.null(BINDING_SITES_MAP[[assay_x]])) {
    data_x <- data_x[!(Pos_real %in% BINDING_SITES_MAP[[assay_x]])]
  }
  if (!is.null(BINDING_SITES_MAP[[assay_y]])) {
    data_y <- data_y[!(Pos_real %in% BINDING_SITES_MAP[[assay_y]])]
  }
  
  data_x_clean <- data_x[, .(mt, Pos_real, ddG, ddG_std)]
  setnames(data_x_clean, "ddG", paste0("ddG_", assay_x))
  setnames(data_x_clean, "ddG_std", paste0("std_", assay_x))
  
  data_y_clean <- data_y[, .(mt, Pos_real, ddG, ddG_std)]
  setnames(data_y_clean, "ddG", paste0("ddG_", assay_y))
  setnames(data_y_clean, "ddG_std", paste0("std_", assay_y))
  
  merged_data <- merge(data_x_clean, data_y_clean, by = c("mt", "Pos_real"))
  
  pvalue_threshold <- function(av, se, threshold){
    zscore <- (abs(av) - threshold) / se
    2 * pnorm(abs(zscore), lower.tail = FALSE)
  }
  
  merged_data[, p_x := pvalue_threshold(get(paste0("ddG_", assay_x)),
                                        get(paste0("std_", assay_x)),
                                        threshold_x)]
  
  merged_data[, p_y := pvalue_threshold(get(paste0("ddG_", assay_y)),
                                        get(paste0("std_", assay_y)),
                                        threshold_y)]
  
  merged_data[, p_adj_x := p.adjust(p_x, "BH")]
  merged_data[, p_adj_y := p.adjust(p_y, "BH")]
  
  list(
    data = merged_data,
    threshold_x = threshold_x,
    threshold_y = threshold_y
  )
}

# =========================================================
# 打印突变列表的辅助函数
# =========================================================
print_mutation_list <- function(mutations, title) {
  cat("\n", title, "\n", sep="")
  cat("突变数:", length(mutations), "\n")
  if(length(mutations) > 0) {
    for(i in seq(1, length(mutations), by=10)) {
      end_idx <- min(i+9, length(mutations))
      cat(paste(mutations[i:end_idx], collapse=", "), "\n")
    }
  } else {
    cat("None\n")
  }
}

# =========================================================
# 通用分析函数（添加NBP标记）
# =========================================================
analyze_protein_pair <- function(protein_x, protein_y,
                                 input_file_x, input_file_y,
                                 anno_file,
                                 verbose = TRUE) {
  
  prepared <- prepare_merged_data_with_FDR(input_file_x, input_file_y, 
                                           protein_x, protein_y, 
                                           anno_file)
  merged_data <- prepared$data
  threshold_x <- prepared$threshold_x
  threshold_y <- prepared$threshold_y
  
  # FDR显著性判断
  merged_data[, pass_FDR_x := p_adj_x < 0.05]
  merged_data[, pass_FDR_y := p_adj_y < 0.05]
  
  # 第一步：基于方向分类（只看阈值）
  merged_data[, direction_class := classify_by_direction(
    get(paste0("ddG_", protein_x)), 
    get(paste0("ddG_", protein_y)), 
    threshold_x, threshold_y
  )]
  
  # 第二步：基于FDR重新评估
  merged_data[, final_classification := reclassify_by_FDR(
    direction_class, pass_FDR_x, pass_FDR_y
  )]
  
  # 设置因子水平
  merged_data[, final_classification := factor(final_classification, levels = LEGEND_ORDER)]
  
  # 标记是否为NBP位点（仅对anticorrelated突变）
  merged_data[, is_NBP := FALSE]
  merged_data[final_classification %in% c("Promoting in X / Disrupting in Y", 
                                          "Disrupting in X / Promoting in Y"), 
              is_NBP := Pos_real %in% NBP_RESIDUES]
  
  if (verbose) {
    cat("\n")
    cat(rep("=", 60), sep="", collapse="")
    cat("\n", protein_x, "vs", protein_y, "分类结果\n")
    cat(rep("=", 60), sep="", collapse="")
    cat("\n")
    
    cat("\n阈值 (", protein_x, "): ", threshold_x, " kcal/mol", sep="")
    cat("\n阈值 (", protein_y, "): ", threshold_y, " kcal/mol", sep="")
    cat("\nFDR阈值: 0.05\n")
    
    # 提取各类突变
    both_promoting <- merged_data[final_classification == "Both promoting", mt]
    both_disrupting <- merged_data[final_classification == "Both disrupting", mt]
    promo_x_disrupt_y <- merged_data[final_classification == "Promoting in X / Disrupting in Y", mt]
    disrupt_x_promo_y <- merged_data[final_classification == "Disrupting in X / Promoting in Y", mt]
    only_x <- merged_data[final_classification == "Allosteric only in X", mt]
    only_y <- merged_data[final_classification == "Allosteric only in Y", mt]
    
    # 统计anticorrelated中NBP位点数量
    anticorrelated_nbp <- merged_data[final_classification %in% c("Promoting in X / Disrupting in Y", 
                                                                  "Disrupting in X / Promoting in Y") & is_NBP == TRUE, .N]
    
    print_mutation_list(both_promoting, paste0("\n1. CORRELATED - Both promoting (", protein_x, " & ", protein_y, "):"))
    print_mutation_list(both_disrupting, paste0("\n2. CORRELATED - Both disrupting (", protein_x, " & ", protein_y, "):"))
    print_mutation_list(promo_x_disrupt_y, paste0("\n3. ANTICORRELATED - Promoting in ", protein_x, " / Disrupting in ", protein_y, ":"))
    print_mutation_list(disrupt_x_promo_y, paste0("\n4. ANTICORRELATED - Disrupting in ", protein_x, " / Promoting in ", protein_y, ":"))
    print_mutation_list(only_x, paste0("\n5. INDEPENDENT - Allosteric only in ", protein_x, ":"))
    print_mutation_list(only_y, paste0("\n6. INDEPENDENT - Allosteric only in ", protein_y, ":"))
    
    cat("\n", rep("=", 60), sep="", collapse="")
    cat("\n统计汇总:\n")
    cat("  总突变数:", nrow(merged_data), "\n")
    cat("  显著突变数:", sum(merged_data$final_classification != "Not significant (FDR >= 0.05)"), "\n")
    cat("  Correlated (促进):", length(both_promoting), "\n")
    cat("  Correlated (破坏):", length(both_disrupting), "\n")
    cat("  Anticorrelated (", protein_x, "促进/", protein_y, "破坏):", length(promo_x_disrupt_y), "\n")
    cat("  Anticorrelated (", protein_x, "破坏/", protein_y, "促进):", length(disrupt_x_promo_y), "\n")
    cat("    - 其中NBP位点:", anticorrelated_nbp, "\n")
    cat("  Independent (仅", protein_x, "):", length(only_x), "\n")
    cat("  Independent (仅", protein_y, "):", length(only_y), "\n")
    cat("  Not significant (FDR >= 0.05):", sum(merged_data$final_classification == "Not significant (FDR >= 0.05)"), "\n")
    cat(rep("=", 60), sep="", collapse="")
    cat("\n")
    
    # 打印分类对比（调试用）
    cat("\n分类对比（方向 vs 最终）:\n")
    comparison <- merged_data[, .N, by = .(direction_class, final_classification)]
    print(comparison)
    
    mutations_list <- list(
      both_promoting = both_promoting,
      both_disrupting = both_disrupting,
      promo_x_disrupt_y = promo_x_disrupt_y,
      disrupt_x_promo_y = disrupt_x_promo_y,
      only_x = only_x,
      only_y = only_y
    )
  } else {
    mutations_list <- NULL
  }
  
  return(list(
    data = merged_data,
    thresholds = c(threshold_x, threshold_y),
    names = c(protein_x, protein_y),
    mutations = mutations_list
  ))
}

# =========================================================
# 通用绘图函数（anticorrelated中NBP用三角形）
# =========================================================
plot_protein_pair <- function(analysis_result, 
                              point_size = 2.5, 
                              alpha = 0.7, 
                              base_size = 15,
                              xlim = c(-1.5, 3),
                              ylim = c(-1.5, 3)) {
  
  merged_data <- analysis_result$data
  protein_x <- analysis_result$names[1]
  protein_y <- analysis_result$names[2]
  threshold_x <- analysis_result$thresholds[1]
  threshold_y <- analysis_result$thresholds[2]
  
  # 计算相关系数
  ddG_x_col <- paste0("ddG_", protein_x)
  ddG_y_col <- paste0("ddG_", protein_y)
  
  cor_test <- cor.test(merged_data[[ddG_x_col]], merged_data[[ddG_y_col]])
  r_value <- round(cor_test$estimate, 3)
  p_value <- cor_test$p.value
  
  sig_label <- ifelse(p_value < 0.001, "***",
                      ifelse(p_value < 0.01, "**",
                             ifelse(p_value < 0.05, "*", "")))
  
  # 创建绘图数据框
  merged_data$final_classification <- factor(merged_data$final_classification, 
                                             levels = LEGEND_ORDER)
  
  # 分离anticorrelated中的NBP和非NBP
  anticorrelated_nbp <- merged_data[final_classification %in% c("Promoting in X / Disrupting in Y", 
                                                                "Disrupting in X / Promoting in Y") & is_NBP == TRUE]
  anticorrelated_non_nbp <- merged_data[final_classification %in% c("Promoting in X / Disrupting in Y", 
                                                                    "Disrupting in X / Promoting in Y") & is_NBP == FALSE]
  
  # 其他显著突变（Both promoting, Both disrupting, Allosteric only）
  other_significant <- merged_data[final_classification %in% c("Both promoting", "Both disrupting",
                                                               "Allosteric only in X", "Allosteric only in Y")]
  
  # 非显著突变
  not_significant <- merged_data[final_classification == "Not significant (FDR >= 0.05)"]
  
  p <- ggplot() +
    theme_classic(base_size = base_size) +
    
    # 阈值线
    geom_vline(xintercept = c(-threshold_x, threshold_x),
               linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_hline(yintercept = c(-threshold_y, threshold_y),
               linetype = "dashed", color = "grey50", linewidth = 0.5) +
    
    # 非显著突变（灰色，较低透明度）- 与参考代码一致
    geom_point(data = not_significant,
               aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"),
               size = point_size, alpha = alpha * 0.3, stroke = 0.3) +
    
    # 其他显著突变（正常透明度）- 与参考代码一致
    geom_point(data = other_significant,
               aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"),
               size = point_size, alpha = alpha, stroke = 0.3) +
    
    # Anticorrelated 非NBP位点（正常透明度，圆形）
    geom_point(data = anticorrelated_non_nbp,
               aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"),
               size = point_size, alpha = alpha, shape = 16, stroke = 0.3) +
    
    # Anticorrelated NBP位点（正常透明度，三角形）
    geom_point(data = anticorrelated_nbp,
               aes_string(x = ddG_x_col, y = ddG_y_col, color = "final_classification"),
               size = point_size + 0.5, alpha = alpha, shape = 17, stroke = 0.3) +
    
    scale_color_manual(values = COLOR_MAP, breaks = LEGEND_ORDER, drop = FALSE) +
    
    # 添加相关系数文本
    annotate("text", x = -Inf, y = Inf, 
             label = paste0("R = ", r_value, sig_label),
             hjust = -0.2, vjust = 1.5, size = base_size/3) +
    
    xlab(paste0("Binding ΔΔG (", protein_x, ") (kcal/mol)")) +
    ylab(paste0("Binding ΔΔG (", protein_y, ") (kcal/mol)")) +
    
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.text = element_text(size = base_size - 2),
      legend.title = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      legend.margin = margin(t = 5, b = 5),
      axis.text = element_text(size = base_size - 2),
      axis.title = element_text(size = base_size),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # 添加这行：x轴文本旋转90度
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(10, 10, 10, 10)
    )+
    
    coord_cartesian(xlim = xlim, ylim = ylim) +
    
    guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 3)))
  
  return(p)
}

# =========================================================
# 使用示例
# =========================================================

# 设置文件路径
base_dir <- "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/"
anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"
output_dir <- "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/"

# 确保输出目录存在
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =========================================================
# 分析 RAF1 vs K13
# =========================================================
cat("\n========== Analyzing RAF1 vs K13 ==========\n")
analysis_RAF1_K13 <- analyze_protein_pair(
  protein_x = "RAF1",
  protein_y = "K13",
  input_file_x = file.path(base_dir, "weights_Binding_RAF.txt"),
  input_file_y = file.path(base_dir, "weights_Binding_K13.txt"),
  anno_file = anno_file,
  verbose = TRUE
)

plot_RAF1_K13 <- plot_protein_pair(analysis_RAF1_K13, 
                                   point_size = 2.5, 
                                   alpha = 0.7, 
                                   base_size = 15)

print(plot_RAF1_K13)

ggsave(file.path(output_dir, "RAF1_vs_K13_scatter_plot_NBP_triangle.pdf"), 
       plot = plot_RAF1_K13, 
       width = 4.5, 
       height = 6,
       device = cairo_pdf)

# =========================================================
# 分析 RAF1 vs K55
# =========================================================
cat("\n========== Analyzing RAF1 vs K55 ==========\n")
analysis_RAF1_K55 <- analyze_protein_pair(
  protein_x = "RAF1",
  protein_y = "K55",
  input_file_x = file.path(base_dir, "weights_Binding_RAF.txt"),
  input_file_y = file.path(base_dir, "weights_Binding_K55.txt"),  
  anno_file = anno_file,
  verbose = TRUE
)

plot_RAF1_K55 <- plot_protein_pair(analysis_RAF1_K55, 
                                   point_size = 2.5, 
                                   alpha = 0.7, 
                                   base_size = 15)

print(plot_RAF1_K55)

ggsave(file.path(output_dir, "RAF1_vs_K55_scatter_plot_NBP_triangle.pdf"), 
       plot = plot_RAF1_K55, 
       width = 4.5, 
       height = 6,
       device = cairo_pdf)

# =========================================================
# 分析 RAF1 vs K27
# =========================================================
cat("\n========== Analyzing RAF1 vs K27 ==========\n")
analysis_RAF1_K27 <- analyze_protein_pair(
  protein_x = "RAF1",
  protein_y = "K27",
  input_file_x = file.path(base_dir, "weights_Binding_RAF.txt"),
  input_file_y = file.path(base_dir, "weights_Binding_K27.txt"),  
  anno_file = anno_file,
  verbose = TRUE
)

plot_RAF1_K27 <- plot_protein_pair(analysis_RAF1_K27, 
                                   point_size = 2.5, 
                                   alpha = 0.7, 
                                   base_size = 15)

print(plot_RAF1_K27)

ggsave(file.path(output_dir, "RAF1_vs_K27_scatter_plot_NBP_triangle.pdf"), 
       plot = plot_RAF1_K27, 
       width = 4.5, 
       height = 6,
       device = cairo_pdf)

# =========================================================
# 分析 K27 vs K13
# =========================================================
cat("\n========== Analyzing K27 vs K13 ==========\n")
analysis_K27_K13 <- analyze_protein_pair(
  protein_x = "K27",
  protein_y = "K13",
  input_file_x = file.path(base_dir, "weights_Binding_K27.txt"),
  input_file_y = file.path(base_dir, "weights_Binding_K13.txt"),
  anno_file = anno_file,
  verbose = TRUE
)

plot_K27_K13 <- plot_protein_pair(analysis_K27_K13, 
                                  point_size = 2.5, 
                                  alpha = 0.7, 
                                  base_size = 15)

print(plot_K27_K13)

ggsave(file.path(output_dir, "K27_vs_K13_scatter_plot_NBP_triangle.pdf"), 
       plot = plot_K27_K13, 
       width = 4.5, 
       height = 6,
       device = cairo_pdf)

# =========================================================
# 分析 K13 vs K19
# =========================================================
cat("\n========== Analyzing K13 vs K19 ==========\n")
analysis_K13_K19 <- analyze_protein_pair(
  protein_x = "K13",
  protein_y = "K19",
  input_file_x = file.path(base_dir, "weights_Binding_K13.txt"),
  input_file_y = file.path(base_dir, "weights_Binding_K19.txt"),
  anno_file = anno_file,
  verbose = TRUE
)

plot_K13_K19 <- plot_protein_pair(analysis_K13_K19, 
                                  point_size = 2.5, 
                                  alpha = 0.7, 
                                  base_size = 15)

print(plot_K13_K19)

ggsave(file.path(output_dir, "K13_vs_K19_scatter_plot_NBP_triangle.pdf"), 
       plot = plot_K13_K19, 
       width = 4.5, 
       height = 6,
       device = cairo_pdf)

