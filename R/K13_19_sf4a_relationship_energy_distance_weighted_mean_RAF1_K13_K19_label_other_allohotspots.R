# ===============================
# Step 0: Load required libraries
# ===============================
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(data.table)
library(dplyr)
library(krasddpcams)
library(usethis)

# ===============================
# Step 1: Define statistical function
# ===============================
print_allosteric_statistics <- function(data_plot, allosteric_list, assays, reg_threshold) {
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("变构位点统计报告\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  # 总体统计
  cat("【总体统计】\n")
  cat(sprintf("回归阈值 (reg_threshold): %.3f kcal/mol\n\n", reg_threshold))
  
  # 按assay统计
  for (assayi in assays) {
    cat(sprintf("【Assay: %s】\n", assayi))
    cat(rep("-", 40), "\n", sep = "")
    
    # 获取当前assay的数据
    data_current <- data_plot[assay == assayi, ]
    
    # 统计各类型位点数量
    binding_sites <- allosteric_list[[assayi]][["Binding interface site"]]
    allosteric_gtp_sites <- allosteric_list[[assayi]][["Allosteric GTP pocket site"]]
    other_gtp_sites <- allosteric_list[[assayi]][["Other GTP pocket site"]]
    major_allosteric_sites <- allosteric_list[[assayi]][["Major allosteric site"]]
    
    # 打印数量统计
    cat("\n位点数量统计：\n")
    cat(sprintf("  • Binding interface sites (结合界面位点): %d\n", length(binding_sites)))
    cat(sprintf("  • Allosteric GTP pocket sites (变构GTP口袋位点): %d\n", length(allosteric_gtp_sites)))
    cat(sprintf("  • Other GTP pocket sites (其他GTP口袋位点): %d\n", length(other_gtp_sites)))
    cat(sprintf("  • Major allosteric sites (主要变构位点): %d\n", length(major_allosteric_sites)))
    
    # 打印具体的残基位置
    if (length(allosteric_gtp_sites) > 0) {
      cat("\n变构GTP口袋位点具体残基：\n")
      cat(sprintf("  残基位置: %s\n", paste(sort(allosteric_gtp_sites), collapse = ", ")))
      
      # 打印这些位点的详细信息
      cat("\n  详细信息 (mean ddG, distance, count):\n")
      for (pos in sort(allosteric_gtp_sites)) {
        pos_data <- data_current[Pos_real == pos, ]
        if (nrow(pos_data) > 0) {
          cat(sprintf("    Pos %d: mean = %.3f, distance = %.2f Å, count = %.1f\n", 
                      pos, pos_data$mean[1], pos_data$distance_bp[1], pos_data$count[1]))
        }
      }
    }
    
    if (length(major_allosteric_sites) > 0) {
      cat("\n主要变构位点具体残基：\n")
      cat(sprintf("  残基位置: %s\n", paste(sort(major_allosteric_sites), collapse = ", ")))
      
      # 打印这些位点的详细信息
      cat("\n  详细信息 (mean ddG, distance, count):\n")
      for (pos in sort(major_allosteric_sites)) {
        pos_data <- data_current[Pos_real == pos, ]
        if (nrow(pos_data) > 0) {
          cat(sprintf("    Pos %d: mean = %.3f, distance = %.2f Å, count = %.1f\n", 
                      pos, pos_data$mean[1], pos_data$distance_bp[1], pos_data$count[1]))
        }
      }
    }
    
    if (length(binding_sites) > 0) {
      cat("\n结合界面位点具体残基：\n")
      cat(sprintf("  残基位置: %s\n", paste(sort(binding_sites), collapse = ", ")))
    }
    
    if (length(other_gtp_sites) > 0) {
      cat("\n其他GTP口袋位点具体残基：\n")
      cat(sprintf("  残基位置: %s\n", paste(sort(other_gtp_sites), collapse = ", ")))
    }
    
    cat("\n")
  }
  
  # 跨assay比较
  cat(rep("=", 80), "\n", sep = "")
  cat("【跨Assay比较】\n")
  cat(rep("-", 40), "\n", sep = "")
  
  # 找出所有assay中共同的和特有的变构位点
  all_allosteric_sites <- list()
  for (assayi in assays) {
    all_allosteric_sites[[assayi]] <- allosteric_list[[assayi]][["Allosteric GTP pocket site"]]
  }
  
  # 交集
  common_sites <- Reduce(intersect, all_allosteric_sites)
  if (length(common_sites) > 0) {
    cat("\n所有assay共有的变构GTP口袋位点：\n")
    cat(sprintf("  残基位置: %s\n", paste(sort(common_sites), collapse = ", ")))
  } else {
    cat("\n所有assay共有的变构GTP口袋位点：无\n")
  }
  
  # 并集
  all_sites <- unique(unlist(all_allosteric_sites))
  cat(sprintf("\n所有assay中出现的变构GTP口袋位点总数: %d\n", length(all_sites)))
  if (length(all_sites) > 0) {
    cat(sprintf("  完整列表: %s\n", paste(sort(all_sites), collapse = ", ")))
  }
  
  # 每个assay特有的位点
  for (i in seq_along(assays)) {
    assayi <- assays[i]
    other_assays <- assays[-i]
    unique_sites <- setdiff(all_allosteric_sites[[assayi]], 
                            unlist(all_allosteric_sites[other_assays]))
    if (length(unique_sites) > 0) {
      cat(sprintf("\n%s特有的变构GTP口袋位点：\n", assayi))
      cat(sprintf("  残基位置: %s\n", paste(sort(unique_sites), collapse = ", ")))
    } else {
      cat(sprintf("\n%s特有的变构GTP口袋位点：无\n", assayi))
    }
  }
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("统计报告结束\n")
  cat(rep("=", 80), "\n\n", sep = "")
}

# ===============================
# Step 2: Main plotting function with cross-protein labeling
# ===============================
plot_weighted_mean_ddG_distance_with_cross <- function(ddG_files, assays, anno_file, 
                                                       x_intercept = 5, 
                                                       output_file = NULL,
                                                       base_font_size = 12,
                                                       point_size = 0.5,
                                                       text_repel_size = 3) {
  
  rect_input <- data.frame(
    xstart = c(3, 15, 38, 51, 67, 77, 87, 109, 127, 139, 148),
    xend = c(9, 24, 44, 57, 73, 84, 104, 115, 136, 143, 166),
    col = c("b1", "a1", "b2", "b3", "a2", "b4", "a3", "b5", "a4", "b6", "a5")
  )
  
  # ===============================
  # Read annotation file
  # ===============================
  anno <- fread(anno_file)
  anno[, Pos_real := Pos]
  
  # ===============================
  # Calculate weighted mean ddG for each assay
  # ===============================
  weighted_mean_ddG <- list()
  for (i in seq_along(assays)) {
    weighted_mean_ddG[[assays[i]]] <- krasddpcams__get_weighted_mean_abs_ddG_mutcount(
      ddG = ddG_files[i],
      assay_sele = assays[i]
    )
  }
  
  # ===============================
  # Combine data for all assays
  # ===============================
  data_plot <- data.table()
  for (assayi in assays) {
    data_plot_assayi <- merge(weighted_mean_ddG[[assayi]], anno, by = "Pos_real", all = TRUE)
    
    data_plot_assayi[, binding_type := "allosteric site"]
    data_plot_assayi[get(paste0("scHAmin_ligand_", assayi)) < x_intercept, 
                     binding_type := "binding site"]
    
    data_plot_assayi[, binding_type_gtp_included := binding_type]
    # Use RAF1 values for GTP binding site determination
    data_plot_assayi[GXPMG_scHAmin_ligand_RAF1 < x_intercept, 
                     binding_type_gtp_included := "GTP binding site"]
    
    data_plot <- rbind(data_plot, data_plot_assayi)
  }
  
  # ===============================
  # Calculate regression threshold using ALL binding sites from ALL assays
  # ===============================
  reg_threshold <- data_plot[binding_type == "binding site", 
                             sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / 
                               sum(1 / .SD[[2]]^2, na.rm = TRUE), 
                             .SDcols = c("mean", "sigma")]
  
  print(paste("Regression threshold (based on all binding sites):", reg_threshold))
  
  # ===============================
  # Classify site types
  # ===============================
  data_plot[, site_type := "Reminder"]
  data_plot[binding_type_gtp_included == "binding site", 
            site_type := "Binding interface site"]
  data_plot[binding_type_gtp_included == "GTP binding site", 
            site_type := "Other GTP pocket site"]
  data_plot[binding_type_gtp_included == "GTP binding site" & 
              mean > reg_threshold & 
              binding_type != "binding site" & 
              count > 9.5, 
            site_type := "Allosteric GTP pocket site"]
  data_plot[binding_type_gtp_included == "allosteric site" & 
              mean > reg_threshold & 
              count > 9.5, 
            site_type := "Major allosteric site"]
  
  # ===============================
  # Annotate secondary structure
  # ===============================
  data_plot[, colors_type := "others"]
  rects_dt <- as.data.table(rect_input)
  
  # Beta strands
  for (b in c("b1", "b2", "b3", "b4", "b5", "b6")) {
    data_plot[Pos_real >= rects_dt[col == b, xstart] & 
                Pos_real <= rects_dt[col == b, xend], 
              colors_type := b]
  }
  
  # Alpha helices
  for (a in c("a1", "a2", "a3", "a4", "a5")) {
    data_plot[Pos_real >= rects_dt[col == a, xstart] & 
                Pos_real <= rects_dt[col == a, xend], 
              colors_type := a]
  }
  
  # Define shapes
  data_plot[, shape := "others"]
  data_plot[colors_type %chin% c("b1", "b2", "b3", "b4", "b5", "b6"), 
            shape := "beta strand"]
  data_plot[colors_type %chin% c("a1", "a2", "a3", "a4", "a5"), 
            shape := "alpha helix"]
  
  # ===============================
  # Filter and set factors
  # ===============================
  data_plot <- data_plot[Pos_real > 1 & count > 9.5, ]
  data_plot <- within(data_plot,
                      site_type <- factor(site_type,
                                          levels = c("Binding interface site",
                                                     "Allosteric GTP pocket site",
                                                     "Other GTP pocket site",
                                                     "Major allosteric site",
                                                     "Reminder")))
  data_plot <- within(data_plot,
                      assay <- factor(assay,
                                      levels = assays))
  
  # ===============================
  # Calculate distances and create allosteric lists
  # ===============================
  allosteric_list <- list()
  for (assayi in assays) {
    data_plot[assay == assayi, 
              distance_bp := get(paste0("scHAmin_ligand_", assayi))]
    allosteric_list[[assayi]] <- list()
    allosteric_list[[assayi]][["Binding interface site"]] <- 
      data_plot[binding_type == "binding site" & assay == assayi, Pos_real]
    allosteric_list[[assayi]][["Allosteric GTP pocket site"]] <- 
      data_plot[site_type == "Allosteric GTP pocket site" & assay == assayi, Pos_real]
    allosteric_list[[assayi]][["Other GTP pocket site"]] <- 
      data_plot[site_type == "Other GTP pocket site" & assay == assayi, Pos_real]
    allosteric_list[[assayi]][["Major allosteric site"]] <- 
      data_plot[site_type == "Major allosteric site" & assay == assayi, Pos_real]
  }
  
  # ===============================
  # Identify cross-protein hotspots (from other proteins, not self)
  # ===============================
  
  protein_names <- assays
  
  # For each protein, identify positions that are allosteric hotspots in OTHER proteins
  # but are NOT already classified as allosteric hotspots in the current protein
  for (current_protein in protein_names) {
    # Get current protein's own allosteric hotspots
    self_allo <- unique(c(
      allosteric_list[[current_protein]][["Allosteric GTP pocket site"]],
      allosteric_list[[current_protein]][["Major allosteric site"]]
    ))
    
    # Get allosteric hotspots from other proteins
    other_proteins <- setdiff(protein_names, current_protein)
    
    # Collect allosteric hotspot positions from other proteins
    other_hotspots <- data.table()
    for (other in other_proteins) {
      other_allo <- unique(c(
        allosteric_list[[other]][["Allosteric GTP pocket site"]],
        allosteric_list[[other]][["Major allosteric site"]]
      ))
      
      if(length(other_allo) > 0) {
        other_hotspots <- rbind(other_hotspots, 
                                data.table(Pos_real = other_allo, 
                                           source_protein = other))
      }
    }
    
    # Remove duplicates and exclude positions that are already self allosteric hotspots
    other_hotspots <- unique(other_hotspots)
    other_hotspots <- other_hotspots[!(Pos_real %in% self_allo), ]
    
    # Mark these positions in data_plot for the current protein
    data_plot[assay == current_protein, 
              cross_hotspot_source := NA_character_]
    
    # Add a column for label color
    data_plot[assay == current_protein, 
              label_color := NA_character_]
    
    for(i in 1:nrow(other_hotspots)) {
      pos <- other_hotspots$Pos_real[i]
      source <- other_hotspots$source_protein[i]
      
      data_plot[assay == current_protein & Pos_real == pos, 
                cross_hotspot_source := source]
      
      # Set label color based on source
      if(source == "K13") {
        data_plot[assay == current_protein & Pos_real == pos, 
                  label_color := "#F1DD10"]
      } else if(source == "K19") {
        data_plot[assay == current_protein & Pos_real == pos, 
                  label_color := "#C68EFD"]
      } else if(source == "RAF1") {
        data_plot[assay == current_protein & Pos_real == pos, 
                  label_color := "#09B636"]
      }
    }
  }
  
  # ===============================
  # Create a unified color mapping for points
  # ===============================
  
  # Assign point colors based on site_type and cross_hotspot_source
  data_plot[, point_color := as.character(site_type)]
  
  # Override colors for cross-hotspots (from other proteins)
  data_plot[!is.na(cross_hotspot_source), 
            point_color := paste0("Cross_", cross_hotspot_source)]
  
  # Define color values
  color_values <- c(
    "Binding interface site" = "#F4270C",      # Red - binding interface
    "Allosteric GTP pocket site" = "#1B38A6",  # Blue - self allosteric GTP pocket
    "Major allosteric site" = "#F4AD0C",       # Yellow - self major allosteric
    "Other GTP pocket site" = "#75C2F6",       # Light blue - other GTP pocket
    "Reminder" = "gray",                       # Gray - others
    "Cross_K13" = "#F1DD10",                   # Red - K13 hotspot (in other proteins)
    "Cross_K19" = "#C68EFD",                   # Blue - K19 hotspot (in other proteins)
    "Cross_RAF1" = "#09B636"                   # Green - RAF1 hotspot (in other proteins)
  )
  
  # Define labels for legend
  color_labels <- c(
    "Binding interface site" = "Binding interface",
    "Allosteric GTP pocket site" = "Allosteric GTP pocket (self)",
    "Major allosteric site" = "Major allosteric (self)",
    "Other GTP pocket site" = "Other GTP pocket",
    "Reminder" = "Others",
    "Cross_K13" = "K13 allosteric hotspot",
    "Cross_K19" = "K19 allosteric hotspot",
    "Cross_RAF1" = "RAF1 allosteric hotspot"
  )
  
  # ===============================
  # Print allosteric statistics
  # ===============================
  print_allosteric_statistics(data_plot, allosteric_list, assays, reg_threshold)
  
  # Print cross-hotspot information
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("【交叉热点标注信息（来自其他蛋白，且非自身变构热点）】\n")
  cat(rep("-", 40), "\n", sep = "")
  
  for (current_protein in protein_names) {
    cross_data <- data_plot[assay == current_protein & !is.na(cross_hotspot_source), ]
    if(nrow(cross_data) > 0) {
      cat(sprintf("\n在%s中标注的其他蛋白热点（非自身变构热点）:\n", current_protein))
      for(i in 1:nrow(cross_data)) {
        cat(sprintf("  残基 %d (来自%s): mean ddG = %.3f, distance = %.2f Å, 点颜色 = %s\n",
                    cross_data$Pos_real[i], 
                    cross_data$cross_hotspot_source[i],
                    cross_data$mean[i],
                    cross_data$distance_bp[i],
                    cross_data$cross_hotspot_source[i]))
      }
    } else {
      cat(sprintf("\n在%s中没有需要额外标注的其他蛋白热点\n", current_protein))
    }
  }
  cat(rep("=", 80), "\n\n", sep = "")
  
  # ===============================
  # Create plot with cross-protein labels and colored points
  # ===============================
  
  p <- ggplot2::ggplot() +
    # Points with unified color mapping
    ggplot2::geom_point(data = data_plot,
                        mapping = aes(x = distance_bp,
                                      y = mean,
                                      color = point_color,
                                      shape = as.factor(shape)),
                        size = point_size) +
    # Error bars
    ggplot2::geom_pointrange(data = data_plot,
                             aes(x = distance_bp,
                                 y = mean,
                                 color = point_color,
                                 ymin = mean - sigma,
                                 ymax = mean + sigma,
                                 shape = as.factor(shape)),
                             size = point_size) +
    # Reference lines
    ggplot2::geom_hline(yintercept = reg_threshold, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = x_intercept, linetype = 2, linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.3) +
    # Text labels for major allosteric sites (self, original logic)
    ggrepel::geom_text_repel(data = data_plot[site_type == "Major allosteric site", ],
                             aes(x = distance_bp,
                                 y = mean,
                                 label = Pos_real),
                             nudge_y = 0.05,
                             color = "#F4AD0C",
                             size = text_repel_size,
                             fontface = "bold") +
    # Text labels for allosteric GTP pocket sites (self, original logic)
    ggrepel::geom_text_repel(data = data_plot[site_type == "Allosteric GTP pocket site", ],
                             aes(x = distance_bp,
                                 y = mean,
                                 label = Pos_real),
                             nudge_y = 0.05,
                             color = "#1B38A6",
                             size = text_repel_size,
                             fontface = "bold") +
    # Text labels for cross-protein hotspots (from other proteins, not self)
    # Using specific colors for each cross hotspot
    ggrepel::geom_text_repel(data = data_plot[!is.na(cross_hotspot_source) & cross_hotspot_source == "K13", ],
                             aes(x = distance_bp,
                                 y = mean,
                                 label = paste0(Pos_real, " (K13)")),
                             nudge_y = 0.08,
                             size = text_repel_size * 0.8,
                             fontface = "bold",
                             color = "#F1DD10") +
    ggrepel::geom_text_repel(data = data_plot[!is.na(cross_hotspot_source) & cross_hotspot_source == "K19", ],
                             aes(x = distance_bp,
                                 y = mean,
                                 label = paste0(Pos_real, " (K19)")),
                             nudge_y = 0.08,
                             size = text_repel_size * 0.8,
                             fontface = "bold",
                             color = "#C68EFD") +
    ggrepel::geom_text_repel(data = data_plot[!is.na(cross_hotspot_source) & cross_hotspot_source == "RAF1", ],
                             aes(x = distance_bp,
                                 y = mean,
                                 label = paste0(Pos_real, " (RAF1)")),
                             nudge_y = 0.08,
                             size = text_repel_size * 0.8,
                             fontface = "bold",
                             color = "#09B636") +
    # Axis labels
    ggplot2::xlab(expression(paste("Distance to binding partner (" * ring(A) * ")"))) +
    ggplot2::ylab("Weighted mean |ddG| (kcal/mol)") +
    ggplot2::labs(color = "Site Type", shape = "Secondary Structure") +
    # Facet
    ggplot2::facet_wrap(~assay, ncol = 3) +
    # Color scale for points (only one scale_color_manual now)
    ggplot2::scale_color_manual(values = color_values,
                                labels = color_labels,
                                breaks = names(color_values),
                                drop = FALSE) +
    # Shape scale
    ggplot2::scale_shape_manual(values = c("beta strand" = 15, 
                                           "alpha helix" = 16, 
                                           "others" = 17),
                                drop = FALSE) +
    # Theme
    ggplot2::theme_classic(base_size = base_font_size) +
    ggplot2::theme(
      axis.text.x = element_text(size = base_font_size * 0.8, vjust = 0.5, hjust = 0.5),
      axis.text.y = element_text(size = base_font_size * 0.8, vjust = 0.5, hjust = 0.5),
      text = element_text(size = base_font_size),
      legend.position = "right",
      strip.text.x = element_text(size = base_font_size),
      strip.background = element_rect(colour = "white", fill = "white"),
      panel.spacing = unit(0.2, "mm"),
      legend.text = element_text(size = base_font_size * 0.6),
      plot.margin = margin(0, 1, 0, 1, "mm"),
      legend.margin = margin(0, 0, 0, -2, "mm"),
      legend.spacing.y = unit(0, 'mm'),
      legend.key.height = unit(4, "mm")
    )
  
  # ===============================
  # Save plot if output file is specified
  # ===============================
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, 
                    plot = p,
                    device = cairo_pdf, 
                    height = 6, 
                    width = 16, 
                    dpi = 300)
    message(paste("Plot saved to:", output_file))
  }
  
  return(list(plot = p, data = data_plot, allosteric_list = allosteric_list, 
              threshold = reg_threshold, cross_hotspots = data_plot[!is.na(cross_hotspot_source), ]))
}

# ===============================
# Step 3: Run the analysis
# ===============================

# Define input files
ddG_files <- c(
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt",
  "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt"
)

assays <- c("K13", "K19", "RAF1")

anno_file <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv"

# Run the function with cross-protein labeling
result <- plot_weighted_mean_ddG_distance_with_cross(
  ddG_files = ddG_files,
  assays = assays,
  anno_file = anno_file,
  x_intercept = 5,
  output_file = "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260406_start version/figure3/K13_K19_RAF1_with_cross_labels_colored_text.pdf",
  base_font_size = 15,
  point_size = 0.8,
  text_repel_size = 5
)

# Display the plot
print(result$plot)

# Access results
print(paste("Threshold value:", result$threshold))

# View cross-hotspot data (only non-self hotspots)
cat("\n\n额外标注的交叉热点（来自其他蛋白，非自身变构热点）:\n")
print(result$cross_hotspots[, .(assay, Pos_real, cross_hotspot_source, mean, distance_bp, site_type)])

# Optional: Save cross-hotspot data to CSV
fwrite(result$cross_hotspots[, .(assay, Pos_real, cross_hotspot_source, mean, sigma, distance_bp, site_type, count)], 
       "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260406_start version/figure3/cross_hotspots_data.csv")
