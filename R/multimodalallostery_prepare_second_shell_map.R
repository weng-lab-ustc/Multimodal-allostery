#' Prepare Second Shell Map.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param contact_shell Value supplied for `contact_shell`.
#' @param assays Value supplied for `assays`.
#' @param shell_value Value supplied for `shell_value`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_second_shell_map <- function(contact_shell, assays, shell_value = 2) {
    result <- setNames(vector("list", length(assays)), assays)
    for (assay in assays) {
        shell_column <- paste0(assay, "_contact_shell")
        result[[assay]] <- if (shell_column %in% names(contact_shell)) {
            contact_shell[get(shell_column) == shell_value, Pos_real]
        }
        else {
            integer()
        }
    }
    result
}

