library(data.table)

wt_aa <- "TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"

# NBP positions (nucleotide-binding pocket)
NBP_positions <- c(12,13,14,15,16,17,18,28,29,30,32,34,35,57,60,61,116,117,119,120,145,146,147)

# Beta sheets positions
beta_sheets <- c(3:9, 38:44, 51:57, 77:84, 109:115, 139:143)

# Function to calculate odds ratio using Fisher's exact test
calculate_odds_ratio <- function(case_in_region, control_in_region, 
                                 case_out_region, control_out_region) {
  contingency_matrix <- matrix(c(case_in_region, control_in_region,
                                 case_out_region, control_out_region), 
                               nrow = 2, byrow = TRUE)
  
  fisher_result <- fisher.test(contingency_matrix)
  odds_ratio <- fisher_result$estimate
  p_value <- fisher_result$p.value
  
  return(list(
    odds_ratio = odds_ratio,
    p_value = p_value,
    matrix = contingency_matrix
  ))
}

# Function to process single assay
process_single_assay <- function(input, assay_sele, anno, contact_shell_file = NULL) {
  
  # Read and process ddG data
  ddG <- fread(input)
  ddG[, Pos_real := Pos_ref + 1]
  ddG[id != "WT", `:=`(wt_codon, substr(id, 1, 1))]
  ddG[id != "WT", `:=`(mt_codon, substr(id, nchar(id), nchar(id)))]
  
  # Create heatmap tool for all possible mutations
  aa_list <- strsplit("GAVLMIFYWKRHDESTCNQP", "")[[1]]
  heatmap_tool <- data.table(wt_codon = rep(strsplit(wt_aa, "")[[1]], each = 20),
                             Pos_real = rep(2:188, each = 20),
                             mt_codon = rep(aa_list, times = length(strsplit(wt_aa, "")[[1]])))
  
  # Merge with ddG data
  ddG <- merge(ddG, heatmap_tool, by = c("Pos_real", "wt_codon", "mt_codon"), all = TRUE)
  ddG[, Pos := Pos_real]
  
  # Calculate weighted mean ddG
  output <- ddG[Pos_real > 1, .(
    mean = sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) / sum(1 / .SD[[2]]^2, na.rm = TRUE)
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  output_sigma <- ddG[Pos_real > 1, .(
    sigma = sqrt(1 / sum(1 / .SD[[2]]^2, na.rm = TRUE))
  ), .SDcols = c("mean_kcal/mol", "std_kcal/mol"), by = "Pos_real"]
  
  weighted_mean_ddG <- merge(output, output_sigma, by = "Pos_real")
  weighted_mean_ddG[, Pos := Pos_real]
  
  # Merge with annotation data
  anno <- fread(anno)
  data_plot <- merge(weighted_mean_ddG, anno, by = "Pos", all = TRUE)
  
  # =========================================================
  # Define binding site types (consistent with Manhattan plot)
  # =========================================================
  
  data_plot[, binding_type := "Reminder"]
  
  # Orthosteric interface
  data_plot[get(paste0("scHAmin_ligand_", assay_sele)) < 5,
            binding_type := "binding site"]
  
  # Include GTP binding interface (assay-specific; same as Manhattan)
  data_plot[, binding_type_gtp_included := binding_type]
  
  data_plot[get(paste0("GXPMG_scHAmin_ligand_", assay_sele)) < 5,
            binding_type_gtp_included := "GTP binding site"]
  
  # =========================================================
  # Calculate regulatory threshold
  # =========================================================
  
  reg_threshold <- data_plot[
    binding_type == "binding site",
    sum(abs(.SD[[1]]) / .SD[[2]]^2, na.rm = TRUE) /
      sum(1 / .SD[[2]]^2, na.rm = TRUE),
    .SDcols = c("mean", "sigma")
  ]
  
  # =========================================================
  # Define site types (consistent with Manhattan plot)
  # =========================================================
  
  data_plot[, site_type := "Reminder"]
  
  data_plot[
    binding_type_gtp_included == "binding site",
    site_type := "Binding interface site"
  ]
  
  data_plot[
    binding_type_gtp_included == "GTP binding site",
    site_type := "GTP binding interface site"
  ]
  
  # =========================================================
  # Merge mutation data with site type
  # =========================================================
  
  data_plot_mutation1 <- merge(
    ddG,
    data_plot[, .(Pos, site_type)],
    by = "Pos",
    all.x = TRUE
  )
  
  data_plot_mutation <- data_plot_mutation1[
    Pos > 1 & !is.na(id)
  ]
  
  data_plot_mutation[, mutation_type := "Reminder"]
  
  # =========================================================
  # Identify allosteric mutations
  # =========================================================
  
  data_plot_mutation[, allosteric_mutation := p.adjust(
    krasddpcams__pvalue(abs(mean) - reg_threshold, std),
    method = "BH"
  ) < 0.05 & (abs(mean) - reg_threshold) > 0]
  
  # =========================================================
  # Direction classification
  # =========================================================
  
  data_plot_mutation[, direction :=
                       ifelse(`mean_kcal/mol` > 0,
                              "inhibit",
                              "stabilize")]
  
  # =========================================================
  # Mutation type classification
  # (same framework as Manhattan plot)
  # =========================================================
  
  data_plot_mutation[
    Pos %in% data_plot[
      site_type == "Binding interface site", Pos
    ] &
      allosteric_mutation == TRUE,
    mutation_type := "Orthosteric site huge differences"
  ]
  
  data_plot_mutation[
    Pos %in% data_plot[
      site_type == "Binding interface site", Pos
    ] &
      allosteric_mutation == FALSE,
    mutation_type := "Orthosteric site small differences"
  ]
  
  data_plot_mutation[
    Pos %in% data_plot[
      site_type == "GTP binding interface site", Pos
    ] &
      allosteric_mutation == TRUE,
    mutation_type := "GTP binding allosteric mutation"
  ]
  
  data_plot_mutation[
    Pos %in% data_plot[
      site_type == "GTP binding interface site", Pos
    ] &
      allosteric_mutation == FALSE,
    mutation_type := "GTP binding other mutation"
  ]
  
  data_plot_mutation[
    !site_type %in% c(
      "GTP binding interface site",
      "Binding interface site"
    ) &
      allosteric_mutation == TRUE,
    mutation_type := "Allosteric mutation"
  ]
  
  data_plot_mutation[
    !site_type %in% c(
      "GTP binding interface site",
      "Binding interface site"
    ) &
      allosteric_mutation == FALSE,
    mutation_type := "Other mutation"
  ]
  
  # =========================================================
  # Define allosteric category
  # =========================================================
  
  data_plot_mutation[, is_allosteric := FALSE]
  
  data_plot_mutation[
    allosteric_mutation == TRUE,
    is_allosteric := TRUE
  ]
  
  # =========================================================
  # Remove orthosteric interface mutations only
  # (keep GTP-binding mutations)
  # =========================================================
  
  mutations_analysis <- data_plot_mutation[
    site_type != "Binding interface site"
  ]
  
  cat("\n", paste0(rep("=", 80), collapse = ""), "\n")
  cat("ENRICHMENT ANALYSIS FOR", toupper(assay_sele), "\n")
  cat("使用 RAF1 定义的 GTP binding site\n")
  cat("排除界面位点 (Binding interface site)\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  
  cat("\n总突变数（非界面）:", nrow(mutations_analysis))
  cat("\n变构突变数:", nrow(mutations_analysis[is_allosteric == TRUE]))
  cat("\n非变构突变数:", nrow(mutations_analysis[is_allosteric == FALSE]), "\n")
  
  # ============================================
  # 1. Enrichment in NBP (nucleotide-binding pocket)
  # ============================================
  
  cat("\n--- 1. Enrichment in Nucleotide-binding Pocket (NBP) ---\n")
  cat("   NBP positions:", paste(NBP_positions, collapse = ", "), "\n")
  cat("   Case: Allosteric mutations\n")
  cat("   Control: Non-allosteric mutations\n")
  cat("   Region: NBP (GTP binding interface site positions)\n")
  cat("   Background: All non-interface positions (excluding NBP)\n\n")
  
  # 标记是否在NBP区域
  mutations_analysis[, in_NBP := Pos %in% NBP_positions]
  
  # 统计
  case_in_NBP <- nrow(mutations_analysis[is_allosteric == TRUE & in_NBP == TRUE])
  control_in_NBP <- nrow(mutations_analysis[is_allosteric == FALSE & in_NBP == TRUE])
  case_out_NBP <- nrow(mutations_analysis[is_allosteric == TRUE & in_NBP == FALSE])
  control_out_NBP <- nrow(mutations_analysis[is_allosteric == FALSE & in_NBP == FALSE])
  
  # 1.1 Combined (inhibit + stabilize)
  result_combined <- calculate_odds_ratio(case_in_NBP, control_in_NBP, case_out_NBP, control_out_NBP)
  
  cat("1.1 All allosteric mutations (combined):\n")
  cat("    OR =", round(result_combined$odds_ratio, 3), "\n")
  cat("    p =", format(result_combined$p_value, scientific = TRUE, digits = 3), "\n")
  cat("    Case in NBP:", case_in_NBP, "\n")
  cat("    Control in NBP:", control_in_NBP, "\n")
  cat("    Case out of NBP:", case_out_NBP, "\n")
  cat("    Control out of NBP:", control_out_NBP, "\n")
  cat("    Contingency table:\n")
  print(result_combined$matrix)
  cat("\n")
  
  # 1.2 Inhibit only
  case_in_NBP_inhibit <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "inhibit" & in_NBP == TRUE])
  control_in_NBP_inhibit <- nrow(mutations_analysis[direction == "inhibit" & is_allosteric == FALSE & in_NBP == TRUE])
  case_out_NBP_inhibit <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "inhibit" & in_NBP == FALSE])
  control_out_NBP_inhibit <- nrow(mutations_analysis[direction == "inhibit" & is_allosteric == FALSE & in_NBP == FALSE])
  
  result_inhibit <- calculate_odds_ratio(case_in_NBP_inhibit, control_in_NBP_inhibit, 
                                         case_out_NBP_inhibit, control_out_NBP_inhibit)
  
  cat("1.2 Inhibit only (ΔΔG > 0):\n")
  cat("    OR =", round(result_inhibit$odds_ratio, 3), "\n")
  cat("    p =", format(result_inhibit$p_value, scientific = TRUE, digits = 3), "\n")
  cat("    Case in NBP:", case_in_NBP_inhibit, "\n")
  cat("    Control in NBP:", control_in_NBP_inhibit, "\n")
  cat("    Case out of NBP:", case_out_NBP_inhibit, "\n")
  cat("    Control out of NBP:", control_out_NBP_inhibit, "\n")
  cat("    Contingency table:\n")
  print(result_inhibit$matrix)
  cat("\n")
  
  # 1.3 Stabilize only
  case_in_NBP_stabilize <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "stabilize" & in_NBP == TRUE])
  control_in_NBP_stabilize <- nrow(mutations_analysis[direction == "stabilize" & is_allosteric == FALSE & in_NBP == TRUE])
  case_out_NBP_stabilize <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "stabilize" & in_NBP == FALSE])
  control_out_NBP_stabilize <- nrow(mutations_analysis[direction == "stabilize" & is_allosteric == FALSE & in_NBP == FALSE])
  
  result_stabilize <- calculate_odds_ratio(case_in_NBP_stabilize, control_in_NBP_stabilize,
                                           case_out_NBP_stabilize, control_out_NBP_stabilize)
  
  cat("1.3 Stabilize only (ΔΔG < 0):\n")
  cat("    OR =", round(result_stabilize$odds_ratio, 3), "\n")
  cat("    p =", format(result_stabilize$p_value, scientific = TRUE, digits = 3), "\n")
  cat("    Case in NBP:", case_in_NBP_stabilize, "\n")
  cat("    Control in NBP:", control_in_NBP_stabilize, "\n")
  cat("    Case out of NBP:", case_out_NBP_stabilize, "\n")
  cat("    Control out of NBP:", control_out_NBP_stabilize, "\n")
  cat("    Contingency table:\n")
  print(result_stabilize$matrix)
  cat("\n")
  
  # ============================================
  # 2. Enrichment in second shell
  # ============================================
  
  second_shell_results <- NULL
  
  if (!is.null(contact_shell_file)) {
    contact_shell <- fread(contact_shell_file)
    
    # Get the appropriate contact shell column
    shell_col <- paste0(assay_sele, "_contact_shell")
    
    if (shell_col %in% names(contact_shell)) {
      # Get second shell positions (contact_shell == 2)
      second_shell_positions <- contact_shell[get(shell_col) == 2, Pos_real]
      # Exclude NBP positions from second shell analysis
      #second_shell_positions <- setdiff(second_shell_positions, NBP_positions)
      
      cat("\n--- 2. Enrichment in Second Shell ---\n")
      cat("   Second shell positions (contact_shell = 2):", 
          paste(head(second_shell_positions, 10), collapse = ", "), 
          if(length(second_shell_positions) > 10) "...", "\n")
      cat("   Number of second shell positions:", length(second_shell_positions), "\n")
      cat("   Case: Allosteric mutations\n")
      cat("   Control: Non-allosteric mutations\n")
      cat("   Region: Second shell positions\n")
      cat("Background: All non-interface positions excluding second shell\n\n")
      
      # 标记是否在second shell区域
      mutations_analysis[, in_second_shell := Pos %in% second_shell_positions]
      
      # 2.1 Combined
      case_in_shell <- nrow(mutations_analysis[is_allosteric == TRUE & in_second_shell == TRUE])
      control_in_shell <- nrow(mutations_analysis[is_allosteric == FALSE & in_second_shell == TRUE])
      case_out_shell <- nrow(mutations_analysis[is_allosteric == TRUE & in_second_shell == FALSE])
      control_out_shell <- nrow(mutations_analysis[is_allosteric == FALSE & in_second_shell == FALSE])
      
      if (case_in_shell + control_in_shell + case_out_shell + control_out_shell > 0) {
        result_shell_combined <- calculate_odds_ratio(case_in_shell, control_in_shell,
                                                      case_out_shell, control_out_shell)
        
        cat("2.1 All allosteric mutations (combined):\n")
        cat("    OR =", round(result_shell_combined$odds_ratio, 3), "\n")
        cat("    p =", format(result_shell_combined$p_value, scientific = TRUE, digits = 3), "\n")
        cat("    Case in second shell:", case_in_shell, "\n")
        cat("    Control in second shell:", control_in_shell, "\n")
        cat("    Case out of second shell:", case_out_shell, "\n")
        cat("    Control out of second shell:", control_out_shell, "\n")
        cat("    Contingency table:\n")
        print(result_shell_combined$matrix)
        cat("\n")
        
        # 2.2 Inhibit only
        case_in_shell_inhibit <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "inhibit" & in_second_shell == TRUE])
        control_in_shell_inhibit <- nrow(mutations_analysis[direction == "inhibit" & is_allosteric == FALSE & in_second_shell == TRUE])
        case_out_shell_inhibit <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "inhibit" & in_second_shell == FALSE])
        control_out_shell_inhibit <- nrow(mutations_analysis[direction == "inhibit" & is_allosteric == FALSE & in_second_shell == FALSE])
        
        result_shell_inhibit <- calculate_odds_ratio(case_in_shell_inhibit, control_in_shell_inhibit,
                                                     case_out_shell_inhibit, control_out_shell_inhibit)
        
        cat("2.2 Inhibit only (ΔΔG > 0):\n")
        cat("    OR =", round(result_shell_inhibit$odds_ratio, 3), "\n")
        cat("    p =", format(result_shell_inhibit$p_value, scientific = TRUE, digits = 3), "\n")
        cat("    Case in second shell:", case_in_shell_inhibit, "\n")
        cat("    Control in second shell:", control_in_shell_inhibit, "\n")
        cat("    Case out of second shell:", case_out_shell_inhibit, "\n")
        cat("    Control out of second shell:", control_out_shell_inhibit, "\n")
        cat("    Contingency table:\n")
        print(result_shell_inhibit$matrix)
        cat("\n")
        
        # 2.3 Stabilize only
        case_in_shell_stabilize <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "stabilize" & in_second_shell == TRUE])
        control_in_shell_stabilize <- nrow(mutations_analysis[direction == "stabilize" & is_allosteric == FALSE & in_second_shell == TRUE])
        case_out_shell_stabilize <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "stabilize" & in_second_shell == FALSE])
        control_out_shell_stabilize <- nrow(mutations_analysis[direction == "stabilize" & is_allosteric == FALSE & in_second_shell == FALSE])
        
        result_shell_stabilize <- calculate_odds_ratio(case_in_shell_stabilize, control_in_shell_stabilize,
                                                       case_out_shell_stabilize, control_out_shell_stabilize)
        
        cat("2.3 Stabilize only (ΔΔG < 0):\n")
        cat("    OR =", round(result_shell_stabilize$odds_ratio, 3), "\n")
        cat("    p =", format(result_shell_stabilize$p_value, scientific = TRUE, digits = 3), "\n")
        cat("    Case in second shell:", case_in_shell_stabilize, "\n")
        cat("    Control in second shell:", control_in_shell_stabilize, "\n")
        cat("    Case out of second shell:", case_out_shell_stabilize, "\n")
        cat("    Control out of second shell:", control_out_shell_stabilize, "\n")
        cat("    Contingency table:\n")
        print(result_shell_stabilize$matrix)
        cat("\n")
        
        # Compile results
        second_shell_results <- data.table(
          Assay = assay_sele,
          Category = c("All_allosteric_combined", "All_allosteric_inhibit", "All_allosteric_stabilize"),
          OR = c(result_shell_combined$odds_ratio, result_shell_inhibit$odds_ratio, result_shell_stabilize$odds_ratio),
          P_value = c(result_shell_combined$p_value, result_shell_inhibit$p_value, result_shell_stabilize$p_value),
          Case_in_region = c(case_in_shell, case_in_shell_inhibit, case_in_shell_stabilize),
          Control_in_region = c(control_in_shell, control_in_shell_inhibit, control_in_shell_stabilize),
          Case_out_region = c(case_out_shell, case_out_shell_inhibit, case_out_shell_stabilize),
          Control_out_region = c(control_out_shell, control_out_shell_inhibit, control_out_shell_stabilize)
        )
      } else {
        cat("2. No mutations found in second shell positions\n\n")
      }
    }
  }
  
  # ============================================
  # 3. Enrichment in Beta Sheets
  # ============================================
  
  cat("\n--- 3. Enrichment in Beta Sheets ---\n")
  cat("   Beta sheet positions:", paste(head(beta_sheets, 10), collapse = ", "), 
      if(length(beta_sheets) > 10) "...", "\n")
  cat("   Number of beta sheet positions:", length(beta_sheets), "\n")
  cat("   Case: Allosteric mutations\n")
  cat("   Control: Non-allosteric mutations\n")
  cat("   Region: Beta sheets\n")
  cat("   Background: All non-interface positions excluding beta sheets\n\n")
  
  # 标记是否在beta sheet区域
  mutations_analysis[, in_beta_sheet := Pos %in% beta_sheets]
  
  # 3.1 Combined
  case_in_beta <- nrow(mutations_analysis[is_allosteric == TRUE & in_beta_sheet == TRUE])
  control_in_beta <- nrow(mutations_analysis[is_allosteric == FALSE & in_beta_sheet == TRUE])
  case_out_beta <- nrow(mutations_analysis[is_allosteric == TRUE & in_beta_sheet == FALSE])
  control_out_beta <- nrow(mutations_analysis[is_allosteric == FALSE & in_beta_sheet == FALSE])
  
  if (case_in_beta + control_in_beta + case_out_beta + control_out_beta > 0) {
    result_beta_combined <- calculate_odds_ratio(case_in_beta, control_in_beta,
                                                 case_out_beta, control_out_beta)
    
    cat("3.1 All allosteric mutations (combined):\n")
    cat("    OR =", round(result_beta_combined$odds_ratio, 3), "\n")
    cat("    p =", format(result_beta_combined$p_value, scientific = TRUE, digits = 3), "\n")
    cat("    Case in beta sheets:", case_in_beta, "\n")
    cat("    Control in beta sheets:", control_in_beta, "\n")
    cat("    Case out of beta sheets:", case_out_beta, "\n")
    cat("    Control out of beta sheets:", control_out_beta, "\n")
    cat("    Contingency table:\n")
    print(result_beta_combined$matrix)
    cat("\n")
    
    # 3.2 Inhibit only
    case_in_beta_inhibit <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "inhibit" & in_beta_sheet == TRUE])
    control_in_beta_inhibit <- nrow(mutations_analysis[direction == "inhibit" & is_allosteric == FALSE & in_beta_sheet == TRUE])
    case_out_beta_inhibit <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "inhibit" & in_beta_sheet == FALSE])
    control_out_beta_inhibit <- nrow(mutations_analysis[direction == "inhibit" & is_allosteric == FALSE & in_beta_sheet == FALSE])
    
    result_beta_inhibit <- calculate_odds_ratio(case_in_beta_inhibit, control_in_beta_inhibit,
                                                case_out_beta_inhibit, control_out_beta_inhibit)
    
    cat("3.2 Inhibit only (ΔΔG > 0):\n")
    cat("    OR =", round(result_beta_inhibit$odds_ratio, 3), "\n")
    cat("    p =", format(result_beta_inhibit$p_value, scientific = TRUE, digits = 3), "\n")
    cat("    Case in beta sheets:", case_in_beta_inhibit, "\n")
    cat("    Control in beta sheets:", control_in_beta_inhibit, "\n")
    cat("    Case out of beta sheets:", case_out_beta_inhibit, "\n")
    cat("    Control out of beta sheets:", control_out_beta_inhibit, "\n")
    cat("    Contingency table:\n")
    print(result_beta_inhibit$matrix)
    cat("\n")
    
    # 3.3 Stabilize only
    case_in_beta_stabilize <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "stabilize" & in_beta_sheet == TRUE])
    control_in_beta_stabilize <- nrow(mutations_analysis[direction == "stabilize" & is_allosteric == FALSE & in_beta_sheet == TRUE])
    case_out_beta_stabilize <- nrow(mutations_analysis[is_allosteric == TRUE & direction == "stabilize" & in_beta_sheet == FALSE])
    control_out_beta_stabilize <- nrow(mutations_analysis[direction == "stabilize" & is_allosteric == FALSE & in_beta_sheet == FALSE])
    
    result_beta_stabilize <- calculate_odds_ratio(case_in_beta_stabilize, control_in_beta_stabilize,
                                                  case_out_beta_stabilize, control_out_beta_stabilize)
    
    cat("3.3 Stabilize only (ΔΔG < 0):\n")
    cat("    OR =", round(result_beta_stabilize$odds_ratio, 3), "\n")
    cat("    p =", format(result_beta_stabilize$p_value, scientific = TRUE, digits = 3), "\n")
    cat("    Case in beta sheets:", case_in_beta_stabilize, "\n")
    cat("    Control in beta sheets:", control_in_beta_stabilize, "\n")
    cat("    Case out of beta sheets:", case_out_beta_stabilize, "\n")
    cat("    Control out of beta sheets:", control_out_beta_stabilize, "\n")
    cat("    Contingency table:\n")
    print(result_beta_stabilize$matrix)
    cat("\n")
    
    # Compile beta sheet results
    beta_sheet_results <- data.table(
      Assay = assay_sele,
      Category = c("All_allosteric_combined", "All_allosteric_inhibit", "All_allosteric_stabilize"),
      OR = c(result_beta_combined$odds_ratio, result_beta_inhibit$odds_ratio, result_beta_stabilize$odds_ratio),
      P_value = c(result_beta_combined$p_value, result_beta_inhibit$p_value, result_beta_stabilize$p_value),
      Case_in_region = c(case_in_beta, case_in_beta_inhibit, case_in_beta_stabilize),
      Control_in_region = c(control_in_beta, control_in_beta_inhibit, control_in_beta_stabilize),
      Case_out_region = c(case_out_beta, case_out_beta_inhibit, case_out_beta_stabilize),
      Control_out_region = c(control_out_beta, control_out_beta_inhibit, control_out_beta_stabilize)
    )
  } else {
    cat("3. No mutations found in beta sheet positions\n\n")
    beta_sheet_results <- NULL
  }
  
  # Compile NBP results
  NBP_results <- data.table(
    Assay = assay_sele,
    Category = c("All_allosteric_combined", "All_allosteric_inhibit", "All_allosteric_stabilize"),
    OR = c(result_combined$odds_ratio, result_inhibit$odds_ratio, result_stabilize$odds_ratio),
    P_value = c(result_combined$p_value, result_inhibit$p_value, result_stabilize$p_value),
    Case_in_NBP = c(case_in_NBP, case_in_NBP_inhibit, case_in_NBP_stabilize),
    Control_in_NBP = c(control_in_NBP, control_in_NBP_inhibit, control_in_NBP_stabilize),
    Case_out_NBP = c(case_out_NBP, case_out_NBP_inhibit, case_out_NBP_stabilize),
    Control_out_NBP = c(control_out_NBP, control_out_NBP_inhibit, control_out_NBP_stabilize)
  )
  
  cat(paste0(rep("=", 80), collapse = ""), "\n\n")
  
  return(list(
    NBP_results = NBP_results,
    second_shell_results = second_shell_results,
    beta_sheet_results = beta_sheet_results
  ))
}

# ============================================
# Run analysis
# ============================================

contact_shell_file <- "C:/Users/36146/OneDrive - USTC/DryLab/Data_analysis_scripts/distance_contacts_analysis/5binder_contact_shell2.csv"

# Process K13
result_K13 <- process_single_assay(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K13.txt",
  assay_sele = "K13",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  contact_shell_file = contact_shell_file
)

# Process K19
result_K19 <- process_single_assay(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_K19.txt",
  assay_sele = "K19",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  contact_shell_file = contact_shell_file
)

# Process RAF1
result_RAF1 <- process_single_assay(
  input = "C:/Users/36146/OneDrive - USTC/DryLab/MoCHI_8binders_l2_e6_RA_old_new_merge_at_mochi_20260104_lr_0.025_2048/task_901/weights/weights_Binding_RAF.txt",
  assay_sele = "RAF1",
  anno = "C:/Users/36146/OneDrive - USTC/DryLab/base_information_for_K13_K19_project/anno_final_for_5.csv",
  contact_shell_file = contact_shell_file
)

# Combine and save results
all_NBP_results <- rbind(
  result_K13$NBP_results,
  result_K19$NBP_results,
  result_RAF1$NBP_results
)

all_second_shell_results <- rbind(
  result_K13$second_shell_results,
  result_K19$second_shell_results,
  result_RAF1$second_shell_results
)

all_beta_sheet_results <- rbind(
  result_K13$beta_sheet_results,
  result_K19$beta_sheet_results,
  result_RAF1$beta_sheet_results
)

output_dir <- "C:/Users/36146/OneDrive - USTC/Manuscripts/K13_K19/figures/20260521_start_Updating/result3_figure3_4/"

write.csv(all_NBP_results, 
          file = paste0(output_dir, "allosteric mutations enrichment_NBP.csv"),
          row.names = FALSE)

if (!is.null(all_second_shell_results) && nrow(all_second_shell_results) > 0) {
  write.csv(all_second_shell_results, 
            file = paste0(output_dir, "allosteric mutations enrichment_second_shell.csv"),
            row.names = FALSE)
}

if (!is.null(all_beta_sheet_results) && nrow(all_beta_sheet_results) > 0) {
  write.csv(all_beta_sheet_results, 
            file = paste0(output_dir, "allosteric mutations enrichment_beta_sheets.csv"),
            row.names = FALSE)
}

cat("\nResults saved to:\n")
cat("  -", paste0(output_dir, "enrichment_NBP.csv\n"))
if (!is.null(all_second_shell_results) && nrow(all_second_shell_results) > 0) {
  cat("  -", paste0(output_dir, "enrichment_second_shell.csv\n"))
}
if (!is.null(all_beta_sheet_results) && nrow(all_beta_sheet_results) > 0) {
  cat("  -", paste0(output_dir, "enrichment_beta_sheets.csv\n"))
}


# Print summary
cat("\n\n", paste0(rep("=", 80), collapse = ""), "\n")
cat("SUMMARY - NBP Enrichment\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
print(all_NBP_results)

if (!is.null(all_second_shell_results) && nrow(all_second_shell_results) > 0) {
  cat("\n\n", paste0(rep("=", 80), collapse = ""), "\n")
  cat("SUMMARY - Second Shell Enrichment\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  print(all_second_shell_results)
}

if (!is.null(all_beta_sheet_results) && nrow(all_beta_sheet_results) > 0) {
  cat("\n\n", paste0(rep("=", 80), collapse = ""), "\n")
  cat("SUMMARY - Beta Sheet Enrichment\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  print(all_beta_sheet_results)
}
