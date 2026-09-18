#' analyze enrichment
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_analyze_enrichment()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param assay_name Input or option used by this workflow; see Usage for the exact interface.
#' @param indep_muts Input or option used by this workflow; see Usage for the exact interface.
#' @param contact_shell Input or option used by this workflow; see Usage for the exact interface.
#' @param anno_file Input or option used by this workflow; see Usage for the exact interface.
#' @param NBP_positions Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_analyze_enrichment <- function (assay_name, indep_muts, contact_shell, anno_file, NBP_positions = c(12, 13, 14, 15, 16, 17, 18, 28, 29, 30, 32, 
    34, 35, 57, 60, 61, 116, 117, 119, 120, 145, 146, 147)) 
{
    cat("\n", paste0(rep("=", 70), collapse = ""), "\n")
    cat("Enrichment Analysis for", assay_name, "\n")
    cat(paste0(rep("=", 70), collapse = ""), "\n")
    extract_pos <- function(mut_string) {
        as.numeric(gsub("[A-Z]", "", mut_string))
    }
    indep_muts[, `:=`(pos, extract_pos(mutation))]
    anno <- .ma_read_table(anno_file)
    interface_col <- paste0("scHAmin_ligand_", assay_name)
    interface_positions <- anno[get(interface_col) < 5, Pos]
    all_positions <- 1:188
    all_positions <- all_positions[!all_positions %in% interface_positions]
    case_positions <- indep_muts$pos
    case_positions <- case_positions[!case_positions %in% interface_positions]
    control_positions <- setdiff(all_positions, case_positions)
    cat("\nPosition counts (non-interface positions only):\n")
    cat(sprintf("  Total positions: %d\n", length(all_positions)))
    cat(sprintf("  Case (independent mutation positions): %d\n", length(case_positions)))
    cat(sprintf("  Control (all other positions): %d\n", length(control_positions)))
    in_nbp <- all_positions %in% NBP_positions
    case_in_nbp <- sum(case_positions %in% NBP_positions)
    case_out_nbp <- length(case_positions) - case_in_nbp
    control_in_nbp <- sum(control_positions %in% NBP_positions)
    control_out_nbp <- length(control_positions) - control_in_nbp
    cat("\n--- NBP Enrichment ---\n")
    cat(sprintf("  Case in NBP: %d\n", case_in_nbp))
    cat(sprintf("  Case out of NBP: %d\n", case_out_nbp))
    cat(sprintf("  Control in NBP: %d\n", control_in_nbp))
    cat(sprintf("  Control out of NBP: %d\n", control_out_nbp))
    res_nbp <- multimodalallostery_fisher_test(case_in = case_in_nbp, case_out = case_out_nbp, control_in = control_in_nbp, 
        control_out = control_out_nbp)
    cat(sprintf("\n  Odds Ratio = %.3f\n", res_nbp$OR))
    cat(sprintf("  P-value = %.2e\n", res_nbp$p))
    shell_col <- paste0(assay_name, "_contact_shell")
    if (shell_col %in% names(contact_shell)) {
        second_shell <- contact_shell[get(shell_col) == 2, Pos_real]
        second_shell <- setdiff(second_shell, interface_positions)
        case_in_shell <- sum(case_positions %in% second_shell)
        case_out_shell <- length(case_positions) - case_in_shell
        control_in_shell <- sum(control_positions %in% second_shell)
        control_out_shell <- length(control_positions) - control_in_shell
        cat("\n--- Second Shell Enrichment ---\n")
        cat(sprintf("  Case in shell: %d\n", case_in_shell))
        cat(sprintf("  Case out of shell: %d\n", case_out_shell))
        cat(sprintf("  Control in shell: %d\n", control_in_shell))
        cat(sprintf("  Control out of shell: %d\n", control_out_shell))
        res_shell <- multimodalallostery_fisher_test(case_in = case_in_shell, case_out = case_out_shell, control_in = control_in_shell, 
            control_out = control_out_shell)
        cat(sprintf("\n  Odds Ratio = %.3f\n", res_shell$OR))
        cat(sprintf("  P-value = %.2e\n", res_shell$p))
    }
    else {
        res_shell <- list(OR = NA, p = NA)
    }
    return(data.table::data.table(Assay = assay_name, Region = c("NBP", "SecondShell"), OR = c(res_nbp$OR, res_shell$OR), 
        P_value = c(res_nbp$p, res_shell$p), Case_in = c(case_in_nbp, if (exists("case_in_shell", inherits = FALSE)) case_in_shell else NA), 
        Case_out = c(case_out_nbp, if (exists("case_out_shell", inherits = FALSE)) case_out_shell else NA), Control_in = c(control_in_nbp, 
            if (exists("control_in_shell", inherits = FALSE)) control_in_shell else NA), Control_out = c(control_out_nbp, 
            if (exists("control_out_shell", inherits = FALSE)) control_out_shell else NA), Case_total = c(length(case_positions), 
            length(case_positions)), Control_total = c(length(control_positions), length(control_positions))))
}

