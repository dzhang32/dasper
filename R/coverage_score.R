#' @describeIn coverage_process Score coverage by their abnormality
#'
#' @export
coverage_score <- function(junctions, coverage, score_func = .zscore, ...) {

    ##### Check user input is correct #####

    if (!identical(names(coverage), c("case", "control"))) {
        stop("coverage should have the names 'case' and 'control'")
    }

    if (!identical(names(coverage[["case"]]), c("exon_coverage_start", "exon_coverage_end", "intron_coverage")) |
        !identical(names(coverage[["control"]]), c("exon_coverage_start", "exon_coverage_end", "intron_coverage"))) {
        stop("coverage matrices should be named 'exon_coverage_start', 'exon_coverage_end', 'intron_coverage'")
    }

    ##### Score coverage in relation to controls #####

    print(stringr::str_c(Sys.time(), " - Generating coverage abnormality score..."))

    coverage_scores <- .coverage_score(coverage, score_func, ...)

    ##### Obtain regions of greatest coverage dysruption #####

    print(stringr::str_c(Sys.time(), " - Obtaining regions with greatest coverage dysruption..."))

    coverage_region_scores_max <- .coverage_score_max(coverage_scores)

    ##### Store output #####

    colnames(coverage_region_scores_max[["regions"]]) <- dimnames(junctions)[[2]]
    colnames(coverage_region_scores_max[["scores"]]) <- dimnames(junctions)[[2]]
    assays(junctions)[["coverage_region"]] <- coverage_region_scores_max[["regions"]]
    assays(junctions)[["coverage_score"]] <- coverage_region_scores_max[["scores"]]

    print(stringr::str_c(Sys.time(), " - done!"))

    return(junctions)
}

#' Calculate coverage abnormality score
#'
#' `.coverage_score` calculates uses `score_func` to generate a score that
#' describes how much coverage in cases deviates from the coverage distribution
#' in controls.
#'
#' @inheritParams coverage_score
#' @inheritParams junction_score
#'
#' @return list containing matrices with values for the normalised coverage.
#'
#' @keywords internal
#' @noRd
.coverage_score <- function(coverage, score_func, ...) {
    regions <- c("exon_coverage_start", "exon_coverage_end", "intron_coverage")
    coverage_scores <- vector(mode = "list", length = length(regions))
    names(coverage_scores) <- regions

    for (region in regions) {
        coverage_score_mat <- matrix(
            ncol = ncol(coverage[["case"]][[region]]),
            nrow = nrow(coverage[["case"]][[region]])
        )

        # for every junction/row score coverage
        for (i in seq_along(coverage_score_mat[, 1])) {
            coverage_score_mat[i, ] <-
                score_func(
                    x = coverage[["case"]][[region]][i, ],
                    y = coverage[["control"]][[region]][i, ],
                    ...
                )
        }

        coverage_scores[[region]] <- coverage_score_mat
    }

    names(coverage_scores) <- names(coverage_scores) %>%
        stringr::str_replace("coverage", "coverage_score")

    return(coverage_scores)
}

#' Obtain region of greatest coverage disruption
#'
#' `.coverage_score_max` will, for each junction take the coverage score which
#' has the highest absolute value. In other words, out of the three regions of
#' interest it will only keep the score for the one across which coverage is
#' most disrupted.
#'
#' @inheritParams coverage_score
#' @inheritParams junction_score
#'
#' @return list containing two matrices. One detailing the regions with the
#'   highest scores and the other with the scores themselves.
#'
#' @keywords internal
#' @noRd
.coverage_score_max <- function(coverage_scores) {
    stopifnot(identical(
        names(coverage_scores),
        c(
            "exon_coverage_score_start",
            "exon_coverage_score_end",
            "intron_coverage_score"
        )
    ))

    coverage_scores_max <- matrix(
        ncol = ncol(coverage_scores[[1]]),
        nrow = nrow(coverage_scores[[1]])
    )

    coverage_region_max <- matrix(
        ncol = ncol(coverage_scores[[1]]),
        nrow = nrow(coverage_scores[[1]])
    )

    # loop across samples/cols
    for (i in seq_along(coverage_scores[["exon_coverage_score_start"]][1, ])) {

        ##### Obtain the most dysregulated region #####

        coverage_scores_per_samp <-
            dplyr::tibble(
                exon_coverage_score_start = coverage_scores[["exon_coverage_score_start"]][, i],
                exon_coverage_score_end = coverage_scores[["exon_coverage_score_end"]][, i],
                intron_coverage_score = coverage_scores[["intron_coverage_score"]][, i]
            )

        # find index of the region with highest absolute score
        coverage_region_max[, i] <- coverage_scores_per_samp %>%
            apply(
                MARGIN = 1,
                FUN = function(x) {
                    which.max(abs(x))
                }
            )

        coverage_scores_max[, i] <- coverage_scores_per_samp %>%
            apply(
                MARGIN = 1,
                FUN = function(x) {
                    x[which.max(abs(x))]
                }
            )
    }

    coverage_region_scores_max <- list(
        regions = coverage_region_max,
        scores = coverage_scores_max
    )

    return(coverage_region_scores_max)
}
