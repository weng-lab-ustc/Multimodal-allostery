#' prepare annotation tracks
#'
#' Workflow-only September 2026 implementation of `multimodalallostery_prepare_annotation_tracks()`.
#'
#' @details Scientific logic follows the current `confirm_scriptcode_202609` workflow. Provenance and old/new comparison are recorded in `FUNCTION_INVENTORY_202609.csv`.
#' @param annotation_tracks Input or option used by this workflow; see Usage for the exact interface.
#' @param correlation_data Input or option used by this workflow; see Usage for the exact interface.
#' @return The prepared data, statistics, or plot described by the function name and Usage.
#' @export
multimodalallostery_prepare_annotation_tracks <- function (annotation_tracks, correlation_data) 
{
    result <- data.table::rbindlist(lapply(names(annotation_tracks), function(track_name) {
        data.table::data.table(track = track_name, Pos_real = annotation_tracks[[track_name]])
    }))
    result <- merge(result, correlation_data[, c("Pos_real", "residue_order")], by = "Pos_real", all.x = TRUE)
    result$track <- factor(result$track, levels = rev(names(annotation_tracks)))
    result$residue_order <- factor(result$residue_order, levels = levels(correlation_data$residue_order))
    result
}

