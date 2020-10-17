#' @describeIn junction_process Score patient junctions by their abnormality
#'
#' @export
junction_score <- function(junctions, score_func = .zscore, ...) {

    ##### Check user input is correct #####

    if (is.null(colData(junctions)[["case_control"]])) {
        stop("Junctions colData (sample metadata) must include the column 'case_control'")
    }

    if (!"norm" %in% names(assays(junctions))) {
        stop("Junctions must include the 'norm' assay")
    }

    ##### Obtaining junction direction #####

    print(stringr::str_c(Sys.time(), " - Calculating the direction of change of junctions..."))

    junctions <- .junction_direction(junctions)

    ##### Calculate junction abnormality score #####

    print(stringr::str_c(Sys.time(), " - Generating junction abnormality score..."))

    junctions <- .junction_score(
        junctions = junctions,
        score_func = score_func,
        ...
    )

    print(stringr::str_c(Sys.time(), " - done!"))

    return(junctions)
}

#' Calculate whether junctions are up or down-regulated
#'
#' `.junction_direction` will label whether junction's counts are up-regulated
#' (+1) or down-regulated (-1) with respect to average control counts.
#'
#' @inheritParams junction_score
#' @param ave_func the function to perform on control junctions to obtain their
#'   average. This is currently by default the [mean][base::mean].
#'
#' @return junctions with additional
#'   [assay][SummarizedExperiment::SummarizedExperiment-class] "direction".
#'
#' @keywords internal
#' @noRd
.junction_direction <- function(junctions, ave_func = mean) {
    # drop = FALSE forces return of matrix
    control_count <- assays(junctions)[["norm"]][, colData(junctions)[["case_control"]] == "control",
        drop = FALSE
    ]

    control_average <- apply(control_count,
        MARGIN = 1,
        FUN = function(x) ave_func(x)
    )

    # perform this across all junctions as the controls needed for
    # .junction_score and yet not removed
    direction <- assays(junctions)[["norm"]] - control_average
    direction <- ifelse(direction > 0, 1, -1)

    assays(junctions)[["direction"]] <- direction

    return(junctions)
}

#' Calculate an abnormality score for each junction
#'
#' `.junction_score` will use `score_func` to score each junction based on how
#' much it's counts deviate from controls.
#'
#' @inheritParams junction_score
#'
#' @return junctions with additional
#'   [assay][SummarizedExperiment::SummarizedExperiment-class] 'score'.
#'
#' @keywords internal
#' @noRd
.junction_score <- function(junctions, score_func = .zscore, ...) {

    # drop = FALSE to force return of matrix
    case_count <- assays(junctions)[["norm"]][, colData(junctions)[["case_control"]] == "case", drop = FALSE]
    control_count <- assays(junctions)[["norm"]][, colData(junctions)[["case_control"]] == "control", drop = FALSE]

    case_score <- matrix(
        nrow = nrow(case_count),
        ncol = ncol(case_count)
    )

    # for the current z-score approach, it would faster to
    # vectorise this, however looping/applying to keep the
    # flexibility for more a complex score_func in future
    for (i in seq_along(junctions)) {
        case_score[i, ] <-
            score_func(
                x = case_count[i, ],
                y = control_count[i, ],
                ...
            )
    }

    # add colnames to prevent dimnames warning from
    # SummarizedExperiment::`assays<-`
    colnames(case_score) <- colnames(case_count)

    # subset for only case samples
    # to enable adding score as an assay
    # otherwise score would have a diff number of columns to junction
    junctions <- junctions[, colData(junctions)[["case_control"]] == "case"]
    assays(junctions)[["score"]] <- case_score

    return(junctions)
}
