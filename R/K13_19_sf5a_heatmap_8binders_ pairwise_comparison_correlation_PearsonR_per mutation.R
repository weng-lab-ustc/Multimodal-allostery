# =========================================================
# 8 Binder Correlation Analysis (Simplified with improved plotting)
# =========================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)

# =========================================================
# 0. 定义 binding interface map（每个 assay 的 interface）
# =========================================================

binding_sites_map <- list(
  RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
  RALGDS = c(24, 25, 31, 33, 36, 37, 38, 39, 40, 41, 56, 64, 67),
  PI3KCG = c(3, 21, 24, 25, 33, 36, 37, 38, 39, 40, 41, 63, 64, 70, 73),
  SOS1 = c(1, 22, 24, 25, 26, 27, 31, 33, 36, 37, 38, 39, 41, 42, 43, 44, 45, 50, 56, 59, 64, 65, 66, 67, 70, 149, 153),
  K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74),
  K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71),
  K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138),
  K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
)

# 定义binder的固定顺序
BINDER_ORDER <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55", "K27", "K13", "K19")

# =========================================================
# 1. 简化的数据读取函数（只读取必要的ddG值）
# =========================================================
load_ddG_data <- function(input, assay_sele) {
  ddG <- fread(input)
  ddG[, Pos_real := Pos_ref + 1]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  ddG[, mt := paste0(wt_codon, Pos_real, mt_codon)]
  ddG <- ddG[id != "WT"]
  
  result <- ddG[, .(mt, Pos_real, `mean_kcal/mol`)]
  setnames(result, "mean_kcal/mol", "ddG")
  result[, assay := assay_sele]
  
  return(result)
}

# =========================================================
# 2. 简化的相关性计算函数
# =========================================================
calculate_correlation_simple <- function(input_x, input_y, assay_x, assay_y) {
  
  # 读取数据
  data_x <- load_ddG_data(input_x, assay_x)
  data_y <- load_ddG_data(input_y, assay_y)
  
  # 移除 binding interface 位点
  if (assay_x %in% names(binding_sites_map)) {
    data_x <- data_x[!(Pos_real %in% binding_sites_map[[assay_x]])]
  }
  if (assay_y %in% names(binding_sites_map)) {
    data_y <- data_y[!(Pos_real %in% binding_sites_map[[assay_y]])]
  }
  
  # 合并共享突变
  merged_data <- merge(
    data_x[, .(mt, Pos_real, ddG)],
    data_y[, .(mt, Pos_real, ddG)],
    by = c("mt", "Pos_real"),
    suffixes = c(paste0("_", assay_x), paste0("_", assay_y))
  )
  
  # 计算 Pearson 相关性
  cor_test <- cor.test(merged_data[[paste0("ddG_", assay_x)]],
                       merged_data[[paste0("ddG_", assay_y)]])
  
  return(list(
    r = round(cor_test$estimate, 3),
    p = cor_test$p.value,
    n = nrow(merged_data)
  ))
}

# =========================================================
# 3. 计算所有 assay 对的相关性（返回矩阵形式）
# =========================================================
calculate_correlation_matrix <- function(input_files, binder_order) {
  assays <- names(input_files)
  n <- length(assays)
  
  # 按照指定顺序重新排列 assays
  assays <- binder_order[binder_order %in% assays]
  n <- length(assays)
  
  # 初始化相关性矩阵
  cor_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(cor_matrix) <- assays
  colnames(cor_matrix) <- assays
  
  p_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(p_matrix) <- assays
  colnames(p_matrix) <- assays
  
  # 计算所有 pairwise 相关性
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        cor_matrix[i, j] <- 1
        p_matrix[i, j] <- 0
      } else if (i < j) {
        cat("Calculating:", assays[i], "vs", assays[j], "\n")
        
        cor_res <- calculate_correlation_simple(
          input_x = input_files[[assays[i]]],
          input_y = input_files[[assays[j]]],
          assay_x = assays[i],
          assay_y = assays[j]
        )
        
        cor_matrix[i, j] <- cor_res$r
        cor_matrix[j, i] <- cor_res$r
        
        p_matrix[i, j] <- cor_res$p
        p_matrix[j, i] <- cor_res$p
      }
    }
  }
  
  return(list(cor_matrix = cor_matrix, p_matrix = p_matrix))
}

# =========================================================
# 4. 绘制聚类热图（使用固定顺序，不聚类）
# =========================================================
plot_correlation_heatmap <- function(cor_matrix, p_matrix, 
                                     output_file = NULL,
                                     width = 10, 
                                     height = 8) {
  
  # 定义分组信息（用于行/列注释）
  BI1_active <- c("RAF1", "RALGDS", "PI3KCG", "SOS1", "K55")
  BI1_inactive <- c("K27")  # 特殊标记的K27
  BI2 <- c("K13", "K19")
  
  # 创建注释数据框，为每个assay指定特定颜色
  # 先创建group标签
  annotation_col <- data.frame(
    Group = ifelse(rownames(cor_matrix) %in% BI1_active, "BI1_active", 
                   ifelse(rownames(cor_matrix) %in% BI2, "BI2", 
                          ifelse(rownames(cor_matrix) %in% BI1_inactive, "BI1_inactive", "Other")))
  )
  rownames(annotation_col) <- rownames(cor_matrix)
  
  # 定义颜色（K27使用pink）
  annotation_colors <- list(
    Group = c(BI1_active = "#F4270C",      # 红色
              BI1_inactive = "pink",  # 粉色
              BI2 = "#1B38A6" )      # 蓝色
              #Other = "white" )    # 白色
    )
    
    # 定义热图颜色（0为白色）
    # 使用更精细的颜色梯度，确保0是白色
    heatmap_colors <- colorRampPalette(c("white", "#F4270C"))(100)
    
    # 创建显著性标记矩阵（显示所有值，包括NS）
    sig_stars <- matrix("", nrow = nrow(p_matrix), ncol = ncol(p_matrix))
    for (i in 1:nrow(p_matrix)) {
      for (j in 1:ncol(p_matrix)) {
        if (i != j) {
          if (p_matrix[i, j] < 0.001) {
            sig_stars[i, j] <- "***"
          } else if (p_matrix[i, j] < 0.01) {
            sig_stars[i, j] <- "**"
          } else if (p_matrix[i, j] < 0.05) {
            sig_stars[i, j] <- "*"
          } else {
            sig_stars[i, j] <- "NS"  # 标记不显著
          }
        }
      }
    }
    
    # 设置对角线为空白
    diag(sig_stars) <- ""
    
    # 绘制热图（不进行聚类，使用固定顺序）
    pheatmap(
      mat = cor_matrix,
      color = heatmap_colors,
      cluster_rows = FALSE,  # 不进行行聚类
      cluster_cols = FALSE,  # 不进行列聚类
      border_color = "grey90",  # 设置边框颜色为白色
      display_numbers = sig_stars,
      number_color = "black",
      number_size = 8,
      fontsize_number = 8,
      main = "Pearson Correlation of Mutational Effects (without BI sites)\nPer mutation",
      xlab = "",
      ylab = "",
      annotation_col = annotation_col,
      annotation_row = annotation_col,
      annotation_colors = annotation_colors,
      show_colnames = TRUE,
      show_rownames = TRUE,
      fontsize_row = 10,
      fontsize_col = 10,
      cellwidth = 40,
      cellheight = 40,
      filename = output_file,
      width = width,
      height = height
    )
}

# =========================================================
# 7. 运行分析
# =========================================================

input_files <- list(
  RAF1   = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
  RALGDS = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAL.txt",
  PI3KCG = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_PI3.txt",
  SOS1   = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_SOS.txt",
  K55    = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K55.txt",
  K27    = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K27.txt",
  K13    = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  K19    = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt"
)

# 7.1 计算相关性矩阵（用于热图）- 使用固定顺序
cat("\n========== Calculating correlation matrix for heatmap ==========\n")
cor_results <- calculate_correlation_matrix(input_files, BINDER_ORDER)
cor_matrix <- cor_results$cor_matrix
p_matrix <- cor_results$p_matrix

# 7.2 输出相关性矩阵
cat("\n========== Correlation Matrix ==========\n")
print(cor_matrix)

# 7.3 绘制热图（使用固定顺序，不聚类，0为白色，NS标记）
cat("\n========== Generating heatmap with fixed order ==========\n")
output_heatmap <- "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result4_figure5_2/correlation_heatmap_8binder_fixed_order 5.pdf"
plot_correlation_heatmap(
  cor_matrix = cor_matrix,
  p_matrix = p_matrix,
  output_file = output_heatmap,
  width = 10,
  height = 8
)

cat("\n========== Analysis Complete ==========\n")
cat("Heatmap saved to:", output_heatmap, "\n")
