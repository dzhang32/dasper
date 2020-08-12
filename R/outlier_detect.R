#' Detecting outlier junctions
#'
#' \code{outlier_detect} will use the features in
#' \code{\link[SummarizedExperiment]{assays}} named \code{feature_names} as
#' input into an unsupervised outlier detection algorithm to score each junction
#' based on how outlier-y it looks in relation to other junctions in the
#' patient. The default \code{score} and \code{outlier_scores} can be
#' calculated using the \code{\link{junction_process}} and
#' \code{\link{coverage_process}} respectively.
#'
#' @inheritParams junction_annot
#' @param feature_names names of assays in \code{junctions} that are to be used
#'   as input into the outlier detection model.
#' @param ... additional arguments passed to the outlier detection model
#'   (isolation forest) for setting parameters.
#'
#' @return junctions as a
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'    object with additional assays named \code{outlier_scores}. The lower the
#'   outlier score, the more of an outlier a junction is observed to be.
#'
#' @seealso for more details on the isolation forest model used:
#'   https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.IsolationForest.html
#'
#' @export
outlier_detect <- function(junctions, feature_names = c("score", "coverage_score"), ...) {

    ##### Check user input #####

    if (!"direction" %in% names(assays(junctions))) {
        stop("junctions must contain a 'direction' assay")
    }

    if (!all(feature_names %in% names(assays(junctions)))) {
        feature_missing <- feature_names[!feature_names %in% names(assays(junctions))]
        stop(stringr::str_c(
            "Assays does not contain the following: ",
            stringr::str_c(feature_missing, collapse = ", ")
        ))
    }

    ##### Score junctions with outlier score #####

    outlier_scores <- matrix(
        nrow = dim(junctions)[1],
        ncol = dim(junctions)[2]
    )

    # for each sample
    for (i in seq_len(dim(junctions)[2])) {
        print(stringr::str_c(
            Sys.time(), " - generating outlier scores for sample ",
            i, "/", dim(junctions)[2]
        ))

        # use data.frame instead of tibble as needed to pass into reticulate
        # add index to preserve order
        features <-
            assays(junctions)[c("direction", feature_names)] %>%
            lapply(FUN = function(x) {
                x[, i]
            }) %>%
            as.data.frame() %>%
            dplyr::mutate(index = dplyr::row_number())

        outlier_scores_samp <- vector(mode = "list", length = 2)

        for (j in seq_along(c(1, -1))) {
            features_up_down <- features %>%
                dplyr::filter(direction == c(1, -1)[j])

            features_up_down[["outlier_score"]] <-
                features_up_down %>%
                dplyr::select(one_of(feature_names)) %>%
                .outlier_score(...)

            outlier_scores_samp[[j]] <- features_up_down
        }

        outlier_scores_samp <- do.call(outlier_scores_samp, what = rbind) %>%
            dplyr::arrange(index)

        stopifnot(identical(outlier_scores_samp[["score"]], features[["score"]]))

        outlier_scores[, i] <- outlier_scores_samp[["outlier_score"]]
    }

    colnames(outlier_scores) <- dimnames(junctions)[[2]]
    assays(junctions)[["outlier_score"]] <- outlier_scores

    print(stringr::str_c(Sys.time(), " - done!"))

    return(junctions)
}
