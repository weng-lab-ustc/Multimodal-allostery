# =========================================================
# 位点特异性突变比较分析（包含WT版本）
# =========================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)

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

# binding interface 位点映射
binding_sites_map <- list(
  RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
  K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74),
  K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71),
  K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138),
  K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
)

# 区域特异性位点
NBP_residues <- c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)
Switch_I_residues <- c(25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40)
Switch_II_residues <- c(58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76)

# =========================================================
# 数据读取函数（修改版 - 为每个位点创建WT行）
# =========================================================
load_mutation_data <- function(input, assay_sele) {
  if (is.character(input)) {
    ddG <- fread(input)
  } else {
    ddG <- as.data.table(input)
  }
  
  ddG[, Pos_real := Pos_ref + 1]
  
  # 处理突变数据
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  
  # 提取突变数据
  mut_data <- ddG[id != "WT", .(
    mt = paste0(wt_codon, Pos_real, mt_codon),
    Pos_real,
    wt_codon,
    mt_codon,
    ddG = `mean_kcal/mol`,
    ddG_std = `std_kcal/mol`
  )]
  
  # 获取所有位点
  all_positions <- unique(mut_data$Pos_real)
  
  # 为每个位点创建WT行（ddG=0, ddG_std=0）
  wt_rows <- data.table(
    mt = paste0("WT", all_positions, "WT"),
    Pos_real = all_positions,
    wt_codon = "WT",
    mt_codon = "WT",
    ddG = 0,
    ddG_std = 0
  )
  
  # 合并突变和WT
  result <- rbind(mut_data, wt_rows)
  
  # 添加assay信息
  result[, assay := assay_sele]
  
  # 按位点和突变类型排序
  setorder(result, Pos_real, mt_codon)
  
  return(result)
}

# =========================================================
# 计算阈值（排除WT）
# =========================================================
calculate_threshold <- function(data, assay_sele, anno) {
  anno_data <- fread(anno)
  
  # 只使用突变数据（排除WT）
  data_no_wt <- data[mt_codon != "WT"]
  
  site_ddG <- data_no_wt[, .(
    mean_ddG = sum(abs(ddG) / ddG_std^2, na.rm = TRUE) / sum(1 / ddG_std^2, na.rm = TRUE)
  ), by = Pos_real]
  
  site_sigma <- data_no_wt[, .(
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
# 第二次分类：基于FDR
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
# 准备合并数据和计算双FDR（包含WT）- 修复版
# =========================================================
prepare_merged_data_with_FDR <- function(input_x, input_y, assay_x, assay_y, anno) {
  data_x <- load_mutation_data(input_x, assay_x)
  data_y <- load_mutation_data(input_y, assay_y)
  
  threshold_x <- calculate_threshold(data_x, assay_x, anno)
  threshold_y <- calculate_threshold(data_y, assay_y, anno)
  
  # 数据清洗（保留 mt_codon）
  data_x_clean <- data_x[, .(mt, Pos_real, mt_codon, ddG, ddG_std)]
  setnames(data_x_clean, "ddG", paste0("ddG_", assay_x))
  setnames(data_x_clean, "ddG_std", paste0("std_", assay_x))
  
  data_y_clean <- data_y[, .(mt, Pos_real, mt_codon, ddG, ddG_std)]
  setnames(data_y_clean, "ddG", paste0("ddG_", assay_y))
  setnames(data_y_clean, "ddG_std", paste0("std_", assay_y))
  
  # 合并数据（包含WT）
  merged_data <- merge(data_x_clean, data_y_clean, by = c("mt", "Pos_real", "mt_codon"), all = TRUE)
  
  # 对于缺失的突变，填充NA（但WT应该在两个数据集中都存在）
  # 填充缺失的ddG为0（对于WT来说，两边都是0）
  merged_data[is.na(get(paste0("ddG_", assay_x))), 
              (paste0("ddG_", assay_x)) := 0]
  merged_data[is.na(get(paste0("ddG_", assay_y))), 
              (paste0("ddG_", assay_y)) := 0]
  merged_data[is.na(get(paste0("std_", assay_x))), 
              (paste0("std_", assay_x)) := 0]
  merged_data[is.na(get(paste0("std_", assay_y))), 
              (paste0("std_", assay_y)) := 0]
  
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
# 添加区域标记函数
# =========================================================
add_region_marker <- function(data) {
  # 先全部初始化为Other
  data[, region := "Other"]
  
  # 优先标记NBP（NBP优先级最高）
  data[Pos_real %in% NBP_residues, region := "NBP"]
  
  # 只对不属于NBP的位点判断是否属于Switch区域
  data[region != "NBP" & Pos_real %in% Switch_I_residues, region := "Switch I"]
  data[region != "NBP" & Pos_real %in% Switch_II_residues, region := "Switch II"]
  
  anticorrelated_types <- c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y")
  
  data[, plot_group := as.character(final_classification)]
  
  data[final_classification %in% anticorrelated_types & region == "NBP", 
       plot_group := "Anticorrelated_NBP"]
  data[final_classification %in% anticorrelated_types & region == "Switch I", 
       plot_group := "Anticorrelated_SwitchI"]
  data[final_classification %in% anticorrelated_types & region == "Switch II", 
       plot_group := "Anticorrelated_SwitchII"]
  data[final_classification %in% anticorrelated_types & region == "Other", 
       plot_group := "Anticorrelated_Other"]
  
  return(data)
}

# =========================================================
# 完整分析函数（所有位点，包含WT）- 修复版
# =========================================================
analyze_all_sites <- function(input_x, input_y, assay_x, assay_y, anno) {
  
  prepared <- prepare_merged_data_with_FDR(
    input_x, input_y,
    assay_x, assay_y,
    anno
  )
  
  merged_data <- prepared$data
  threshold_x <- prepared$threshold_x
  threshold_y <- prepared$threshold_y
  
  # =====================================================
  # 确保 mt_codon 列存在
  # =====================================================
  # 从 mt 列提取突变类型
  if (!"mt_codon" %in% colnames(merged_data)) {
    # 如果 mt_codon 不存在，从 mt 列提取
    merged_data[, mt_codon := ifelse(
      grepl("^WT.*WT$", mt),  # 如果 mt 以 WT 开头和结尾
      "WT",
      substr(mt, nchar(mt), nchar(mt))  # 否则取最后一个字符
    )]
  }
  
  # =====================================================
  # FDR判断（WT的p值设为NA，因为std=0）
  # =====================================================
  
  # 先创建列
  merged_data[, `:=`(
    pass_FDR_x = NA,
    pass_FDR_y = NA
  )]
  
  # WT的p值设为NA，不参与FDR计算
  merged_data[mt_codon == "WT", `:=`(
    pass_FDR_x = NA,
    pass_FDR_y = NA,
    p_adj_x = NA,
    p_adj_y = NA
  )]
  
  # 突变数据正常计算
  merged_data[mt_codon != "WT", `:=`(
    pass_FDR_x = p_adj_x < 0.05,
    pass_FDR_y = p_adj_y < 0.05
  )]
  
  # =====================================================
  # 第一步：方向分类（WT分类为"WT"）
  # =====================================================
  
  merged_data[, direction_class := NA_character_]
  merged_data[mt_codon == "WT", direction_class := "WT"]
  
  merged_data[mt_codon != "WT", direction_class := classify_by_direction(
    get(paste0("ddG_", assay_x)), 
    get(paste0("ddG_", assay_y)), 
    threshold_x,
    threshold_y
  )]
  
  # =====================================================
  # 第二步：FDR重新分类（WT保持为"WT"）
  # =====================================================
  
  merged_data[, final_classification := NA_character_]
  merged_data[mt_codon == "WT", final_classification := "WT"]
  
  merged_data[mt_codon != "WT", final_classification := reclassify_by_FDR(
    direction_class,
    pass_FDR_x,
    pass_FDR_y
  )]
  
  # =====================================================
  # 设置factor顺序（添加WT）
  # =====================================================
  
  legend_order_with_WT <- c("WT", LEGEND_ORDER)
  
  merged_data[, final_classification := factor(
    final_classification,
    levels = legend_order_with_WT
  )]
  
  # =====================================================
  # 添加区域标记
  # =====================================================
  
  merged_data <- add_region_marker(merged_data)
  
  # WT的plot_group设为"WT"
  merged_data[mt_codon == "WT", plot_group := "WT"]
  
  # =====================================================
  # 保存阈值
  # =====================================================
  
  threshold_vector <- setNames(
    c(threshold_x, threshold_y),
    c(assay_x, assay_y)
  )
  
  # =====================================================
  # 返回结果
  # =====================================================
  
  return(list(
    data = merged_data,
    thresholds = threshold_vector,
    assays = c(assay_x, assay_y)
  ))
}

# =========================================================
# 修改后的绘图函数（包含WT点）
# =========================================================
plot_site_scatter <- function(full_analysis_result, 
                              target_position,
                              point_size = 4, 
                              alpha = 0.8, 
                              base_size = 14,
                              show_labels = TRUE,
                              label_all = TRUE,
                              xlim = NULL,
                              ylim = NULL,
                              show_WT = TRUE,  # 新增：是否显示WT点
                              WT_point_size = 3,  # WT点大小
                              WT_color = "black",  # WT点颜色
                              WT_shape = 19,  # WT点形状（实心圆）
                              WT_label = "WT"  # WT标签
) {
  
  # 从完整结果中筛选指定位点
  full_data <- full_analysis_result$data
  site_data <- full_data[Pos_real == target_position]
  
  if (nrow(site_data) == 0) {
    stop(paste("Position", target_position, "not found in the data"))
  }
  
  assay_x <- full_analysis_result$assays[1]
  assay_y <- full_analysis_result$assays[2]
  threshold_x <- full_analysis_result$thresholds[assay_x]
  threshold_y <- full_analysis_result$thresholds[assay_y]
  
  # 计算相关系数（包含WT）
  cor_data <- site_data[mt_codon != "WT"]  # 先只用突变数据计算
  cor_test <- cor.test(cor_data[[paste0("ddG_", assay_x)]], 
                       cor_data[[paste0("ddG_", assay_y)]])
  r_value <- round(cor_test$estimate, 3)
  p_value <- cor_test$p.value
  
  # 计算包含WT的相关系数（用于对比）
  cor_data_with_WT <- site_data
  cor_test_with_WT <- cor.test(cor_data_with_WT[[paste0("ddG_", assay_x)]], 
                               cor_data_with_WT[[paste0("ddG_", assay_y)]])
  r_value_with_WT <- round(cor_test_with_WT$estimate, 3)
  p_value_with_WT <- cor_test_with_WT$p.value
  
  # 定义形状映射
  site_data[, plot_shape := "regular"]
  site_data[final_classification == "Not significant (FDR >= 0.05)", plot_shape := "Not significant"]
  site_data[final_classification == "WT", plot_shape := "WT"]
  
  other_types <- c("Both promoting", "Both disrupting", 
                   "Allosteric only in X", "Allosteric only in Y")
  site_data[final_classification %in% other_types, plot_shape := "other_significant"]
  
  # 为anticorrelated突变设置不同形状
  anticorrelated_types <- c("Promoting in X / Disrupting in Y", "Disrupting in X / Promoting in Y")
  site_data[final_classification %in% anticorrelated_types & plot_group == "Anticorrelated_NBP", 
            plot_shape := "NBP"]
  site_data[final_classification %in% anticorrelated_types & plot_group == "Anticorrelated_SwitchI", 
            plot_shape := "SwitchI"]
  site_data[final_classification %in% anticorrelated_types & plot_group == "Anticorrelated_SwitchII", 
            plot_shape := "SwitchII"]
  site_data[final_classification %in% anticorrelated_types & plot_group == "Anticorrelated_Other", 
            plot_shape := "Other_anticorrelated"]
  
  # 形状值（添加WT形状）
  shape_values <- c(
    "regular" = 16,
    "Not significant" = 16,
    "other_significant" = 16,
    "NBP" = 17,
    "SwitchI" = 18,
    "SwitchII" = 15,
    "Other_anticorrelated" = 8,
    "WT" = WT_shape  # WT使用指定形状
  )
  
  # 颜色值（添加WT颜色）
  color_values <- c(
    "Not significant (FDR >= 0.05)" = "grey80",
    "Other (neutral in both)" = "grey80",
    "Both promoting" = "#FFB0A5",           
    "Both disrupting" = "#F4270C",          
    "Promoting in X / Disrupting in Y" = "#F4AD0C",  
    "Disrupting in X / Promoting in Y" = "#F1DD10",  
    "Allosteric only in X" = "#1B38A6",     
    "Allosteric only in Y" = "#75C2F6"      
  )
  
  # =========================================================
  # 确定需要标记的突变
  # =========================================================
  if (label_all) {
    # 标记所有突变（不包括WT）
    data_to_label <- site_data[mt_codon != "WT"]
    cat("\n标记所有突变，共", nrow(data_to_label), "个\n")
  } else {
    # 只标记anticorrelated突变
    anticorrelated_groups <- c("Anticorrelated_NBP", "Anticorrelated_SwitchI", 
                               "Anticorrelated_SwitchII", "Anticorrelated_Other")
    data_to_label <- site_data[plot_group %in% anticorrelated_groups]
    cat("\n只标记anticorrelated突变，共", nrow(data_to_label), "个\n")
  }
  
  # =========================================================
  # 分离WT和突变数据
  # =========================================================
  wt_data <- site_data[mt_codon == "WT"]
  mut_data <- site_data[mt_codon != "WT"]
  
  # 创建图形
  p <- ggplot() +
    
    theme_classic(base_size = base_size) +
    
    # 阈值线
    geom_vline(xintercept = c(-threshold_x, threshold_x),
               linetype = "dashed", color = "grey60", linewidth = 0.8) +
    geom_hline(yintercept = c(-threshold_y, threshold_y),
               linetype = "dashed", color = "grey60", linewidth = 0.8) +
    
    # 非显著突变
    geom_point(data = mut_data[final_classification == "Not significant (FDR >= 0.05)"],
               aes(x = .data[[paste0("ddG_", assay_x)]], 
                   y = .data[[paste0("ddG_", assay_y)]],
                   color = final_classification, 
                   shape = plot_shape),
               size = point_size * 0.7, alpha = alpha * 0.3, stroke = 0.3) +
    
    # 显著突变
    geom_point(data = mut_data[final_classification != "Not significant (FDR >= 0.05)"],
               aes(x = .data[[paste0("ddG_", assay_x)]], 
                   y = .data[[paste0("ddG_", assay_y)]],
                   color = final_classification, 
                   shape = plot_shape),
               size = point_size, alpha = alpha, stroke = 0.8)
  
  # 添加WT点（如果show_WT为TRUE且WT数据存在）
  if (show_WT && nrow(wt_data) > 0) {
    p <- p + geom_point(
      data = wt_data,
      aes(x = .data[[paste0("ddG_", assay_x)]], 
          y = .data[[paste0("ddG_", assay_y)]],
          shape = plot_shape),
      color = WT_color,
      size = WT_point_size,
      stroke = 1.2
    )
    
    # 添加WT标签
    if (show_labels) {
      p <- p + geom_text_repel(
        data = wt_data,
        aes(x = .data[[paste0("ddG_", assay_x)]], 
            y = .data[[paste0("ddG_", assay_y)]],
            label = WT_label),
        color = WT_color,
        size = 4,
        fontface = "bold",
        box.padding = 0.3,
        point.padding = 0.2,
        segment.color = "grey30",
        segment.size = 0.5,
        segment.alpha = 0.8,
        min.segment.length = 0,
        force = 1,
        force_pull = 0.5,
        seed = 123,
        show.legend = FALSE
      )
    }
  }
  
  # 继续添加其他图层
  p <- p +
    scale_color_manual(values = color_values, breaks = names(color_values), drop = FALSE) +
    scale_shape_manual(values = shape_values, breaks = names(shape_values), drop = FALSE) +
    
    # 相关系数标注（显示包含WT和不包含WT的R值）
    annotate("text", x = -Inf, y = Inf, 
             label = paste0("R = ", r_value_with_WT,
                            ifelse(p_value_with_WT < 0.001, "***",
                                   ifelse(p_value_with_WT < 0.01, "**",
                                          ifelse(p_value_with_WT < 0.05, "*", " ns")))),
             hjust = -0.1, vjust = 1.5, size = base_size/3.5) +
    
    labs(
      x = bquote(Binding~Delta*Delta*G~"(" * .(assay_x) * ") (kcal/mol)"),
      y = bquote(Binding~Delta*Delta*G~"(" * .(assay_y) * ") (kcal/mol)"),
      title = paste0(assay_x, " vs ", assay_y, " - Position ", target_position)
      #subtitle = paste0("WT (0,0) point shown in ", WT_color)
    ) +
    
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
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # X轴刻度旋转90度
      axis.title = element_text(size = base_size),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(10, 10, 10, 10),
      plot.title = element_text(hjust = 0.5, size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1, color = "grey40")
    )
  
  # 设置坐标轴范围
  if (!is.null(xlim) && !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim, clip = "off")
  } else {
    x_range <- range(site_data[[paste0("ddG_", assay_x)]], na.rm = TRUE)
    y_range <- range(site_data[[paste0("ddG_", assay_y)]], na.rm = TRUE)
    x_pad <- diff(x_range) * 0.15
    y_pad <- diff(y_range) * 0.15
    p <- p + coord_cartesian(
      xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
      ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
      clip = "off"
    )
  }
  
  # =========================================================
  # 添加突变标签
  # =========================================================
  if (show_labels && nrow(data_to_label) > 0) {
    p <- p + geom_text_repel(
      data = data_to_label,
      aes(
        x = .data[[paste0("ddG_", assay_x)]],
        y = .data[[paste0("ddG_", assay_y)]],
        label = mt,
        color = final_classification
      ),
      size = 3.5,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = "grey50",
      segment.size = 0.3,
      segment.alpha = 0.6,
      min.segment.length = 0,
      max.overlaps = Inf,
      force = 1,
      force_pull = 0.5,
      seed = 123,
      show.legend = FALSE
    )
  }
  
  # 图例设置
  p <- p + guides(
    color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 3)),
    shape = "none"
  )
  
  return(p)
}

# =========================================================
# 执行分析
# =========================================================

anno <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"

# 配置文件路径
input_files <- list(
  RAF1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
  K13  = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt"
)

# 运行完整分析（包含WT）
full_result <- analyze_all_sites(
  input_x = input_files$RAF1,
  input_y = input_files$K13,
  assay_x = "RAF1",
  assay_y = "K13",
  anno = anno
)

# 检查WT是否被正确包含
cat("\n=== 数据验证 ===\n")
cat("总行数:", nrow(full_result$data), "\n")
cat("WT行数:", nrow(full_result$data[mt_codon == "WT"]), "\n")
cat("位点数:", length(unique(full_result$data$Pos_real)), "\n")
#cat("每个位点的WT:", full_result$data[mt_codon == "WT", .N, by = Pos_real], "\n")

# =========================================================
# 绘制位点（包含WT）
# =========================================================

# 绘制位点145
p_pos145 <- plot_site_scatter(
  full_analysis_result = full_result,
  target_position = 145,
  point_size = 5,
  base_size = 20,
  show_labels = TRUE,
  label_all = TRUE,
  xlim = c(-1.5, 3.3),
  ylim = c(-1.5, 3.3),
  show_WT = TRUE,
  WT_point_size = 5,
  WT_color = "black",
  WT_shape = 19,
  WT_label = "WT"
)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result5_figure6/RAF1_vs_K13_scatter_plot_position_145_with_WT.pdf", 
       plot = p_pos145, width = 6.5, height = 8, device = cairo_pdf)

# 绘制位点15
p15 <- plot_site_scatter(
  full_analysis_result = full_result,
  target_position = 15,
  point_size = 5,
  base_size = 20,
  show_labels = TRUE,
  label_all = TRUE,
  xlim = c(-1.5, 3.3),
  ylim = c(-1.5, 3.3),
  show_WT = TRUE,
  WT_point_size = 5,
  WT_color = "black",
  WT_shape = 19,
  WT_label = "WT"
)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result5_figure6/RAF1_vs_K13_scatter_plot_position_15_with_WT.pdf", 
       plot = p15, width = 6.5, height = 8, device = cairo_pdf)

# 绘制位点48
p48 <- plot_site_scatter(
  full_analysis_result = full_result,
  target_position = 48,
  point_size = 5,
  base_size = 20,
  show_labels = TRUE,
  label_all = TRUE,
  xlim = c(-1.5, 3.3),
  ylim = c(-1.5, 3.3),
  show_WT = TRUE,
  WT_point_size = 5,
  WT_color = "black",
  WT_shape = 19,
  WT_label = "WT"
)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result5_figure6/RAF1_vs_K13_scatter_plot_position_48_with_WT.pdf", 
       plot = p48, width = 6.5, height = 8, device = cairo_pdf)

# 绘制位点55
p55 <- plot_site_scatter(
  full_analysis_result = full_result,
  target_position = 55,
  point_size = 5,
  base_size = 20,
  show_labels = TRUE,
  label_all = TRUE,
  xlim = c(-1.5, 3.3),
  ylim = c(-1.5, 3.3),
  show_WT = TRUE,
  WT_point_size = 5,
  WT_color = "black",
  WT_shape = 19,
  WT_label = "WT"
)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result5_figure6/RAF1_vs_K13_scatter_plot_position_55_with_WT.pdf", 
       plot = p55, width = 6.5, height = 8, device = cairo_pdf)

# 绘制位点77
p77 <- plot_site_scatter(
  full_analysis_result = full_result,
  target_position = 77,
  point_size = 5,
  base_size = 20,
  show_labels = TRUE,
  label_all = TRUE,
  xlim = c(-1.5, 3.3),
  ylim = c(-1.5, 3.3),
  show_WT = TRUE,
  WT_point_size = 5,
  WT_color = "black",
  WT_shape = 19,
  WT_label = "WT"
)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result5_figure6/RAF1_vs_K13_scatter_plot_position_77_with_WT.pdf", 
       plot = p77, width = 6.5, height = 8, device = cairo_pdf)

# 绘制位点163
p163 <- plot_site_scatter(
  full_analysis_result = full_result,
  target_position = 163,
  point_size = 5,
  base_size = 20,
  show_labels = TRUE,
  label_all = TRUE,
  xlim = c(-1.5, 3.3),
  ylim = c(-1.5, 3.3),
  show_WT = TRUE,
  WT_point_size = 5,
  WT_color = "black",
  WT_shape = 19,
  WT_label = "WT"
)

ggsave("C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result5_figure6/RAF1_vs_K13_scatter_plot_position_163_with_WT.pdf", 
       plot = p163, width = 6.5, height = 8, device = cairo_pdf)

cat("\n所有图片已保存完成！\n")
cat("绘制的位点：15, 48, 55, 77, 145, 163\n")
cat("每个散点图都包含WT (0,0) 点，用黑色实心圆标记\n")
cat("图中显示两个R值：不含WT和包含WT的相关系数\n")

