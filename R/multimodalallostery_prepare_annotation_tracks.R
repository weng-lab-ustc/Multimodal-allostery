#' Prepare Annotation Tracks.
#'
#' Reusable data preparation logic consolidated from the selected Figure/Panel implementations.
#'
#' @param annotation_tracks Value supplied for `annotation_tracks`.
#' @param correlation_data Value supplied for `correlation_data`.
#'
#' @return A prepared data object used by downstream analysis.
#' @export
multimodalallostery_prepare_annotation_tracks <- function(annotation_tracks, correlation_data) {
    result <- data.table::rbindlist(lapply(names(annotation_tracks), function(track_name) {
        data.table::data.table(track = track_name, Pos_real = annotation_tracks[[track_name]])
    }))
    result <- merge(result, correlation_data[, c("Pos_real", "residue_order")], by = "Pos_real", all.x = TRUE)
    result$track <- factor(result$track, levels = rev(names(annotation_tracks)))
    result$residue_order <- factor(result$residue_order, levels = levels(correlation_data$residue_order))
    result
}

