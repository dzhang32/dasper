#' Wrapper function for processing coverage
#'
#' \code{junction_process} wraps all "coverage_" prefixed functions in
#' \code{dasper}. This is designed to simplify processing of the covearge data
#' for those familiar or uninterested with the intermediates.
#'
#' @inheritParams coverage_norm
#' @inheritParams coverage_score
#'
#' @return
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   object containing junction data with coverage scores.
#'
#' @examples
#'
#' @export
coverage_process <- function(junctions,
    ref,
    unannot_width = 20,
    coverage_paths_case,
    coverage_paths_control,
    coverage_chr_control = NULL,
    load_func = .coverage_load,
    norm_const = 1,
    score_func = .zscore,
    ...) {
    print("# Loading and normalising coverage ---------------------------------------------")

    coverage <- coverage_norm(junctions,
        ref,
        unannot_width = 20,
        coverage_paths_case = coverage_paths_case,
        coverage_paths_control = coverage_paths_control,
        coverage_chr_control = NULL,
        load_func = load_func,
        norm_const = norm_const
    )

    print("# Scoring coverage ---------------------------------------------")

    junctions <- coverage_score(junctions, coverage, score_func, ...)

    return(junctions)
}
