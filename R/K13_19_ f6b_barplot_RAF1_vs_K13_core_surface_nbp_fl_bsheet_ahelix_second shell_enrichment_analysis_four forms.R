library(data.table)
library(ggplot2)
library(dplyr)

# =========================================================
# pair sequence（X axis）
# =========================================================
PAIR_ORDER <- c(
  "RAF1 vs K55",
  "RAF1 vs K27",
  "K55 vs K27",
  "RAF1 vs K13",
  "RAF1 vs K19",
  "K55 vs K13",
  "K55 vs K19",
  "K27 vs K13",
  "K27 vs K19",
  "K13 vs K19"
)

# binding interface residuea
binding_sites_map <- list(
  RAF1 = c(21, 25, 29, 31, 33, 36, 37, 38, 39, 40, 41, 67, 71),
  K55 = c(5, 24, 25, 31, 33, 36, 37, 38, 39, 40, 54, 56, 64, 66, 67, 70, 73, 74),
  K27 = c(21, 24, 25, 27, 31, 33, 36, 38, 39, 40, 41, 43, 52, 54, 67, 70, 71),
  K13 = c(63, 68, 87, 88, 90, 91, 92, 94, 95, 96, 97, 98, 99, 101, 102, 105, 106, 107, 129, 133, 136, 137, 138),
  K19 = c(68, 87, 88, 90, 91, 92, 94, 95, 97, 98, 99, 101, 102, 105, 107, 108, 125, 129, 133, 136, 137)
)



# =========================================================
# function for reading data
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
# First classification: Based on the direction of the effect (considering only the threshold, not the FDR).
# =========================================================
classify_by_direction <- function(ddG_x, ddG_y, threshold_x, threshold_y) {
  sig_x <- abs(ddG_x) > threshold_x
  sig_y <- abs(ddG_y) > threshold_y
  
  disrupt_x <- ddG_x > threshold_x
  promote_x <- ddG_x < -threshold_x
  disrupt_y <- ddG_y > threshold_y
  promote_y <- ddG_y < -threshold_y
  
  result <- rep("neutral", length(ddG_x))
  
  # Correlated/Anticorrelated type (doubly significant)
  result[sig_x & sig_y & promote_x & promote_y] <- "Both promoting"
  result[sig_x & sig_y & disrupt_x & disrupt_y] <- "Both disrupting"
  result[sig_x & sig_y & promote_x & disrupt_y] <- "Promoting in X / Disrupting in Y"
  result[sig_x & sig_y & disrupt_x & promote_y] <- "Disrupting in X / Promoting in Y"
  
  # "Allosteric-only" type (singly significant)
  result[sig_x & !sig_y] <- "Allosteric only in X"
  result[!sig_x & sig_y] <- "Allosteric only in Y"
  
  # Not at all significant
  result[!sig_x & !sig_y] <- "Not significant (FDR >= 0.05)"
  
  return(result)
}

# =========================================================
# Second classification: Re-evaluation of all types based on FDR
# =========================================================
reclassify_by_FDR <- function(direction_class, pass_FDR_x, pass_FDR_y) {
  result <- direction_class
  
  for (i in 1:length(direction_class)) {
    
    # Handling correlated/anticorrelated types (requires dual-FDR significance)
    if (direction_class[i] %in% c("Both promoting", "Both disrupting", 
                                  "Promoting in X / Disrupting in Y", 
                                  "Disrupting in X / Promoting in Y")) {
      if (!(pass_FDR_x[i] & pass_FDR_y[i])) {
        # Downgrade: Check if the single FDR criterion is met.
        if (pass_FDR_x[i] & !pass_FDR_y[i]) {
          result[i] <- "Allosteric only in X"
        } else if (!pass_FDR_x[i] & pass_FDR_y[i]) {
          result[i] <- "Allosteric only in Y"
        } else {
          result[i] <- "Not significant (FDR >= 0.05)"
        }
      }
    }
    
    # Handle the "Allosteric only in X" type (requires significant FDR in X).
    else if (direction_class[i] == "Allosteric only in X") {
      if (!pass_FDR_x[i]) {
        result[i] <- "Not significant (FDR >= 0.05)"
      }
    }
    
    # Handle "Allosteric-only in Y" types (requires significant FDR in Y).
    else if (direction_class[i] == "Allosteric only in Y") {
      if (!pass_FDR_y[i]) {
        result[i] <- "Not significant (FDR >= 0.05)"
      }
    }
  }
  
  return(result)
}

# =========================================================
# Calculate threshold
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
# function: Prepare to merge data and calculate FDR.
# =========================================================
prepare_merged_data_with_FDR <- function(input_x, input_y, assay_x, assay_y, anno) {
  data_x <- load_mutation_data(input_x, assay_x)
  data_y <- load_mutation_data(input_y, assay_y)
  
  threshold_x <- calculate_threshold(data_x, assay_x, anno)
  threshold_y <- calculate_threshold(data_y, assay_y, anno)
  
  if (!is.null(binding_sites_map[[assay_x]])) {
    data_x <- data_x[!(Pos_real %in% binding_sites_map[[assay_x]])]
  }
  if (!is.null(binding_sites_map[[assay_y]])) {
    data_y <- data_y[!(Pos_real %in% binding_sites_map[[assay_y]])]
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
# Main function for two-step classification
# =========================================================
classify_two_step <- function(merged_data, threshold_x, threshold_y, assay_x, assay_y) {
  
  ddG_x_col <- paste0("ddG_", assay_x)
  ddG_y_col <- paste0("ddG_", assay_y)
  
  # Step 1: Classify based on the direction of the effect (consider only the threshold, not the FDR).
  merged_data[, direction_class := classify_by_direction(
    get(ddG_x_col), get(ddG_y_col), threshold_x, threshold_y
  )]
  
  # Calculate FDR flags
  merged_data[, pass_FDR_x := p_adj_x < 0.05]
  merged_data[, pass_FDR_y := p_adj_y < 0.05]
  
  # Step 2: Re-evaluate all types based on FDR.
  merged_data[, final_classification := reclassify_by_FDR(direction_class, pass_FDR_x, pass_FDR_y)]
  
  return(merged_data)
}

# =========================================================
# Calculate the enrichment odds ratio (OR) for each structural region.
# =========================================================
calc_or_original <- function(df, region_residues, cat) {
  df[, in_region := Pos_real %in% region_residues]
  df[, is_cat := category == cat]
  
  a <- sum(df$in_region & df$is_cat)
  b <- sum(!df$in_region & df$is_cat)
  c <- sum(df$in_region & !df$is_cat)
  d <- sum(!df$in_region & !df$is_cat)
  
  tab <- matrix(c(a, b, c, d), nrow = 2)
  
  ft <- fisher.test(tab)
  
  list(OR = unname(ft$estimate), 
       p = ft$p.value, 
       OR_low = ft$conf.int[1], 
       OR_high = ft$conf.int[2])
}

# =========================================================
# Add the second-shell site definition (before `structure_regions`).
# =========================================================

# Read second-shell data
contact_shell_file <- "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/distance_contacts_analysis/5binder_contact_shell2.csv"
contact_shell <- fread(contact_shell_file)

second_shell_map <- list()

for(assay in names(input_files)) {
  shell_col <- paste0(assay, "_contact_shell")
  if(shell_col %in% names(contact_shell)) {
    second_shell_map[[assay]] <- contact_shell[get(shell_col) == 2, Pos_real]
    cat(assay, ":", length(second_shell_map[[assay]]), "second shell positions\n")
  } else {
    second_shell_map[[assay]] <- c()
    cat(assay, ": No contact shell column found\n")
  }
}


run_pair_all_regions <- function(pair_name, input_files, anno, structure_regions, second_shell_map) {
  
  assays <- strsplit(pair_name, " vs ")[[1]]
  x <- assays[1]
  y <- assays[2]
  
  cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
  cat("Analyzing:", pair_name, "\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  
  prepared <- prepare_merged_data_with_FDR(
    input_files[[x]], input_files[[y]], x, y, anno
  )
  
  df <- prepared$data
  threshold_x <- prepared$threshold_x
  threshold_y <- prepared$threshold_y
  
  # Use a two-step classification logic.
  df <- classify_two_step(df, threshold_x, threshold_y, x, y)
  
  # Create a category based on final_classification.
  df[, category := fifelse(
    final_classification %in% c("Both promoting", "Both disrupting"),
    "Correlated",
    fifelse(
      final_classification %in% c("Promoting in X / Disrupting in Y",
                                  "Disrupting in X / Promoting in Y"),
      "Anti-correlated",
      fifelse(
        final_classification == "Allosteric only in X",
        "Allosteric only in X",
        fifelse(
          final_classification == "Allosteric only in Y",
          "Allosteric only in Y",
          "Other"
        )
      )
    )
  )]
  
  # Print statistics
  cat("\n分类统计:\n")
  cat("  Total mutations:", nrow(df), "\n")
  cat("  Correlated:", sum(df$category == "Correlated"), "\n")
  cat("  Anti-correlated:", sum(df$category == "Anti-correlated"), "\n")
  cat("  Allosteric only in X:", sum(df$category == "Allosteric only in X"), "\n")
  cat("  Allosteric only in Y:", sum(df$category == "Allosteric only in Y"), "\n")
  cat("  Other:", sum(df$category == "Other"), "\n")
  
  # Calculate the proportion and odds ratio (OR) for each structural region.
  all_results <- list()
  
  for(region_name in names(structure_regions)) {
    region_residues <- structure_regions[[region_name]]
    
    # Calculate the proportion of each category.
    plot_df <- df[, .(
      frac = mean(Pos_real %in% region_residues),
      n = .N,
      se = sqrt(mean(Pos_real %in% region_residues) * (1 - mean(Pos_real %in% region_residues)) / .N)
    ), by = category]
    plot_df[, region := region_name]
    plot_df[, pair := pair_name]
    
    # Calculate the OR value for each category.
    categories_to_test <- c("Correlated", "Anti-correlated", 
                            "Allosteric only in X", "Allosteric only in Y", "Other")
    or_list <- lapply(categories_to_test, function(cat) {
      res <- calc_or_original(df, region_residues, cat)
      data.table(
        pair = pair_name,
        region = region_name,
        category = cat,
        OR = res$OR,
        OR_low = res$OR_low,
        OR_high = res$OR_high,
        p = res$p
      )
    })
    or_df <- rbindlist(or_list)
    
    all_results[[region_name]] <- list(plot = plot_df, or = or_df)
  }
  
  # =========================================================
  # Add second-shell analysis (using assay-specific second shells).
  # =========================================================
  
  cat("\n--- Second Shell Analysis ---\n")
  
  # 获取X和Y各自的second shell positions
  second_shell_x <- second_shell_map[[x]]
  second_shell_y <- second_shell_map[[y]]
  
  cat("  ", x, "second shell positions:", length(second_shell_x), "\n")
  cat("  ", y, "second shell positions:", length(second_shell_y), "\n")
  
  # =========================================================
  # 对X的second shell进行分析（只对Allosteric only in X）
  # =========================================================
  if(length(second_shell_x) > 0) {
    cat("\n  Analyzing", x, "second shell for 'Allosteric only in X'\n")
    
    # 只保留 Allosteric only in X 和 Other 类别
    df_subset_x <- df[category %in% c("Allosteric only in X", "Other")]
    
    # 计算占比
    plot_df_shell_x <- df_subset_x[, .(
      frac = mean(Pos_real %in% second_shell_x),
      n = .N,
      se = sqrt(mean(Pos_real %in% second_shell_x) * (1 - mean(Pos_real %in% second_shell_x)) / .N)
    ), by = category]
    plot_df_shell_x[, region := paste0("Second Shell (", x, ")")]
    plot_df_shell_x[, pair := pair_name]
    
    # 计算OR值（只对 Allosteric only in X vs Other）
    cat_x <- "Allosteric only in X"
    res_x <- calc_or_original(df_subset_x, second_shell_x, cat_x)
    or_df_shell_x <- data.table(
      pair = pair_name,
      region = paste0("Second Shell (", x, ")"),
      category = cat_x,
      OR = res_x$OR,
      OR_low = res_x$OR_low,
      OR_high = res_x$OR_high,
      p = res_x$p
    )
    
    # 添加Other的OR（为1，作为参考）
    or_df_shell_x_other <- data.table(
      pair = pair_name,
      region = paste0("Second Shell (", x, ")"),
      category = "Other",
      OR = 1,
      OR_low = 1,
      OR_high = 1,
      p = 1
    )
    or_df_shell_x <- rbind(or_df_shell_x, or_df_shell_x_other)
    
    all_results[[paste0("Second Shell (", x, ")")]] <- list(
      plot = plot_df_shell_x, 
      or = or_df_shell_x
    )
  }
  
  # =========================================================
  # 对Y的second shell进行分析（只对Allosteric only in Y）
  # =========================================================
  if(length(second_shell_y) > 0) {
    cat("\n  Analyzing", y, "second shell for 'Allosteric only in Y'\n")
    
    # 只保留 Allosteric only in Y 和 Other 类别
    df_subset_y <- df[category %in% c("Allosteric only in Y", "Other")]
    
    # 计算占比
    plot_df_shell_y <- df_subset_y[, .(
      frac = mean(Pos_real %in% second_shell_y),
      n = .N,
      se = sqrt(mean(Pos_real %in% second_shell_y) * (1 - mean(Pos_real %in% second_shell_y)) / .N)
    ), by = category]
    plot_df_shell_y[, region := paste0("Second Shell (", y, ")")]
    plot_df_shell_y[, pair := pair_name]
    
    # 计算OR值（只对 Allosteric only in Y vs Other）
    cat_y <- "Allosteric only in Y"
    res_y <- calc_or_original(df_subset_y, second_shell_y, cat_y)
    or_df_shell_y <- data.table(
      pair = pair_name,
      region = paste0("Second Shell (", y, ")"),
      category = cat_y,
      OR = res_y$OR,
      OR_low = res_y$OR_low,
      OR_high = res_y$OR_high,
      p = res_y$p
    )
    
    # 添加Other的OR（为1，作为参考）
    or_df_shell_y_other <- data.table(
      pair = pair_name,
      region = paste0("Second Shell (", y, ")"),
      category = "Other",
      OR = 1,
      OR_low = 1,
      OR_high = 1,
      p = 1
    )
    or_df_shell_y <- rbind(or_df_shell_y, or_df_shell_y_other)
    
    all_results[[paste0("Second Shell (", y, ")")]] <- list(
      plot = plot_df_shell_y, 
      or = or_df_shell_y
    )
  }
  
  # 合并所有结果
  combined_plot <- rbindlist(lapply(all_results, `[[`, "plot"))
  combined_or <- rbindlist(lapply(all_results, `[[`, "or"))
  
  cat("\n最终区域列表:", paste(unique(combined_plot$region), collapse=", "), "\n")
  
  list(plot = combined_plot, or = combined_or, full_data = df)
}

# =========================================================
# 定义 input_files 和 anno（在读取 second shell 之前）
# =========================================================
input_files <- list(
  RAF1 = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
  K55  = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K55.txt",
  K27  = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K27.txt",
  K13  = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  K19  = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt"
)

anno <- "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_8.csv"

# =========================================================
# 定义结构区域（在调用 run_pair_all_regions 之前）
# =========================================================

core_residues <- c(4,6,7,8,9,10,11,14,15,16,17,18,19,20,21,22,23,24,40,
                   42,44,46,51,52,53,54,55,56,57,58,68,71,72,75,77,78,79,
                   80,81,82,83,84,89,90,92,93,96,97,99,100,101,103,109,110,
                   111,112,113,114,115,116,118,125,130,133,134,137,139,141,
                   142,143,144,145,146,151,152,155,156,157,158,159,160,162,163)

NBP <- c(12,13,14,15,16,17,18,28,29,30,32,34,35,57,60,61,116,117,119,120,145,146,147)

functional_loop_residues <- c(
  10, 11, 12, 13, 14, 15, 16, 17, 
  25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 
  58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76
)

beta_sheets <- c(3:9, 38:44, 51:57, 77:84, 109:115, 139:143)

surface_residues <- c(1,2,3,5,12,13,25:39,41,43,45,47:50,59:67,
                      69,70,73,74,76,85:88,91,94,95,98,102,104:108,
                      117,119:124,126:129,131,132,135,136,138,140,
                      147:150,153,154,161,164:188)

alpha_helices <- c(15:24, 67:73, 87:104, 127:136, 148:166)

# =========================================================
# 定义 structure_regions
# =========================================================
structure_regions <- list(
  "Core" = core_residues,
  "Surface" = surface_residues,
  "NBP" = NBP,
  "Functional Loop" = functional_loop_residues,
  "Beta Sheets" = beta_sheets,
  "Alpha Helices" = alpha_helices
)

# =========================================================
# 添加Second shell位点定义（在input_files定义之后）
# =========================================================

# 读取Second shell数据
contact_shell_file <- "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/distance_contacts_analysis/5binder_contact_shell2.csv"
contact_shell <- fread(contact_shell_file)

# 为每个assay定义second shell positions (contact_shell == 2)
second_shell_map <- list()

for(assay in names(input_files)) {
  shell_col <- paste0(assay, "_contact_shell")
  if(shell_col %in% names(contact_shell)) {
    second_shell_map[[assay]] <- contact_shell[get(shell_col) == 2, Pos_real]
    cat(assay, ":", length(second_shell_map[[assay]]), "second shell positions\n")
  } else {
    second_shell_map[[assay]] <- c()
    cat(assay, ": No contact shell column found\n")
  }
}

# 打印second_shell_map信息
cat("\n=== Second Shell Map ===\n")
for(assay in names(second_shell_map)) {
  cat(assay, ":", length(second_shell_map[[assay]]), "positions\n")
  if(length(second_shell_map[[assay]]) > 0) {
    cat("  First 10:", paste(head(second_shell_map[[assay]], 10), collapse=", "), "\n")
  }
}

# =========================================================
# 现在可以安全地调用 run_pair_all_regions
# =========================================================
target_pair <- "RAF1 vs K13"

res <- run_pair_all_regions(target_pair, input_files, anno, structure_regions, second_shell_map)

plot_df <- res$plot
or_df <- res$or

# =========================================================
# 创建图：只包含 Allosteric only in X, Allosteric only in Y, Other
# =========================================================

# 设置要显示的类别（只显示Allosteric和Other）
categories_to_plot <- c("Allosteric only in X", "Allosteric only in Y", "Other")

# 颜色映射
color_map <- c(
  "Allosteric only in X" = "#C68EFD", 
  "Allosteric only in Y" = "#09B636",  
  "Other" = "grey80"
)

# 筛选数据
plot_subset <- plot_df[category %in% categories_to_plot]

# 获取所有唯一的区域
all_regions <- unique(plot_subset$region)
cat("\n所有区域:", paste(all_regions, collapse=", "), "\n")

# 设置区域顺序（Second Shell放在最后）
region_levels <- c("Core", "Surface", "NBP", "Functional Loop", "Beta Sheets", "Alpha Helices")
# 添加Second Shell区域
second_shell_regions <- grep("Second Shell", all_regions, value = TRUE)
region_levels <- c(region_levels, sort(second_shell_regions))

plot_subset[, region := factor(region, levels = region_levels)]
plot_subset[, category := factor(category, levels = categories_to_plot)]

# 获取OR值用于标注
or_subset <- or_df[pair == target_pair & category %in% categories_to_plot]
or_subset <- or_subset[region %in% region_levels]
or_subset[, region := factor(region, levels = region_levels)]

# 为每个类别设置x偏移量（3个类别）
or_subset[, x_offset := dplyr::case_when(
  category == "Allosteric only in X" ~ -0.25,
  category == "Allosteric only in Y" ~ 0,
  category == "Other" ~ 0.25
)]

# 设置固定的Y轴位置放置OR标签
fixed_y_position <- 1.2

# 计算Y轴上限
max_y <- max(plot_subset$frac + plot_subset$se, na.rm = TRUE)
if (fixed_y_position > max_y) {
  max_y <- fixed_y_position + 0.1
} else {
  max_y <- max_y * 1.15
}
max_y <- max(max_y, 1.0)

# 创建标签文本
or_subset[, label := paste0("OR = ", round(OR, 2), 
                            ifelse(p < 0.05, "*", ""))]

# 提取assays名称
assays <- strsplit(target_pair, " vs ")[[1]]
x_name <- assays[1]
y_name <- assays[2]

# 创建主图
p <- ggplot(plot_subset, aes(x = region, y = frac, fill = category)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  
  # error bar
  geom_errorbar(aes(ymin = frac - se, ymax = frac + se),
                position = position_dodge(0.9), 
                width = 0.15, 
                size = 0.8,
                color = "#F1DD10",
                alpha = 0.8,
                na.rm = TRUE) +
  
  # error bar端点小横线
  geom_errorbar(aes(ymin = frac - se, ymax = frac - se),
                position = position_dodge(0.9),
                width = 0.3,
                size = 0.8,
                color = "#F1DD10",
                alpha = 0.8,
                na.rm = TRUE) +
  geom_errorbar(aes(ymin = frac + se, ymax = frac + se),
                position = position_dodge(0.9),
                width = 0.3,
                size = 0.8,
                color = "#F1DD10",
                alpha = 0.8,
                na.rm = TRUE) +
  
  # OR标签
  geom_text(data = or_subset, 
            aes(x = as.numeric(region) + x_offset, 
                y = fixed_y_position, 
                label = label,
                color = category),
            size = 3.5,
            angle = 45,
            hjust = 0.5,
            vjust = 0) +
  
  scale_fill_manual(values = color_map,
                    labels = c(paste0("Allosteric only in ", x_name),
                               paste0("Allosteric only in ", y_name),
                               "Other")) +
  scale_color_manual(values = color_map,
                     labels = c(paste0("Allosteric only in ", x_name),
                                paste0("Allosteric only in ", y_name),
                                "Other"),
                     guide = "none") +
  
  scale_y_continuous(limits = c(0, max_y), 
                     expand = expansion(mult = c(0, 0.02)),
                     breaks = seq(0, 1, 0.2)) +
  
  # y = 1 的虚线
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.5, size = 0.8) +
  
  labs(y = "Fraction of mutations in region",
       x = "Structural region") +
  
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 15),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 40, r = 10, b = 10, l = 10)
  )

# 显示图
print(p)

# =========================================================
# 保存图片
# =========================================================
output_dir <- "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result5_figure6/"

ggsave(paste0(output_dir, "barplot_RAF1_vs_K13_allosteric_only_with_second_shell.pdf"),
       plot = p,
       width = 8,
       height = 6,
       device = cairo_pdf)

# =========================================================
# 保存OR值表格
# =========================================================
or_output <- or_df[, .(pair, region, category, OR, OR_low, OR_high, p)]
or_output[, p_signif := ifelse(p < 0.05, "*", "ns")]
or_output[, OR_text := paste0(round(OR, 2), ifelse(p < 0.05, "*", ""))]

fwrite(or_output, 
       paste0(output_dir, "OR_values_RAF1_vs_K13_allosteric_only_with_second_shell.csv"))

# 打印Second Shell统计信息
cat("\n\n========================================\n")
cat("Second Shell Analysis Summary\n")
cat("========================================\n")

for(assay in names(second_shell_map)) {
  cat(assay, ":", length(second_shell_map[[assay]]), "positions\n")
  if(length(second_shell_map[[assay]]) > 0) {
    cat("  Positions:", paste(sort(second_shell_map[[assay]])[1:min(10, length(second_shell_map[[assay]]))], 
                              collapse=", "), 
        if(length(second_shell_map[[assay]]) > 10) "...", "\n")
  }
}

# 打印Second Shell的OR结果
cat("\nSecond Shell OR values:\n")
shell_or <- or_output[grepl("Second Shell", region)]
if(nrow(shell_or) > 0) {
  print(shell_or[, .(region, category, OR_text, p_signif)])
} else {
  cat("  No Second Shell results found.\n")
}
