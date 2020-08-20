#' Wrapper for processing outliers
#'
#' `outlier_process` wraps all "outlier_" prefixed functions in
#' [dasper]. This is designed to simplify processing of the detecting
#' outlier junctions for those familiar or uninterested with the intermediates.
#'
#' @inheritParams junction_annot
#' @inheritParams outlier_detect
#' @inheritParams outlier_aggregate
#'
#' @return `DataFrame` with one row per cluster detailing each cluster's
#'   associated junctions, outlier scores, ranks and genes.
#'
#' @family outlier
#' @export
outlier_process <- function(junctions,
    feature_names = c("score", "coverage_score"),
    samp_id_col = "samp_id",
    ...) {
    print("# Detecting outliers using an isolation forest ------------------------------")

    junctions <- outlier_detect(
        junctions,
        feature_names,
        ...
    )

    print("# Aggregating outlier scores to cluster-level -----------------------------")

    outlier_scores_tidy <- outlier_aggregate(
        junctions,
        samp_id_col
    )

    return(outlier_scores_tidy)
}
