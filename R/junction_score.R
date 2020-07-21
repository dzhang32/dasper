#' Score patient junctions in terms of abnormality in relation to controls
#'
#' \code{junction_score} will use the counts contained within the "norm"
#' \code{\link[SummarizedExperiment]{assay}} to calculate a deviation of each
#' patient junction from the expected distribution of controls counts. The
#' function used to calculate this abnormality score can be user-inputted or
#' left as the default z-score.
#'
#' @inheritParams junction_annot
#'
#' @param score_func function to score junctions by abnromality. By default,
#'   will use a z-score. But this can be switched dependent on the format of the
#'   user's junction data. Function must take as input an x and y argument. x
#'   being the case junctions and y being the controls. This must return a
#'   numeric vector equal to the length of x with each element corresponding to
#'   a abnormality of a junction in each case sample.
#' @param ... additional arguments passed to \code{score_func}.
#'
#' @return junctions as a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}} object filtered
#'   for only "case" samples with an additional \code{assay} containing junction
#'   abnormality scores.
#'
#' @examples
#'
#' junctions <- junction_norm(junctions_example)
#' junctions <- junction_score(junctions)
#' junctions
#' @export
junction_score <- function(junctions, score_func = .zscore, ...) {

    ##### Check user input is correct #####

    if (is.null(colData(junctions)[["case_control"]])) {
        stop("Junctions colData (sample metadata) must include the column 'case_control'")
    }

    if (!"norm" %in% names(assays(junctions))) {
        stop("Junctions must include the 'norm' assay")
    }

    ##### Calculate junction abnormality score #####

    print(stringr::str_c(Sys.time(), " - Generating junction abnormality score..."))

    case_count <- assays(junctions)[["norm"]][, colData(junctions)[["case_control"]] == "case"]
    control_count <- assays(junctions)[["norm"]][, colData(junctions)[["case_control"]] == "control"]

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

    # subset for only case samples
    # to enable adding score as an assay
    # otherwise score would have a diff number of columns to junction
    junctions <- junctions[, colData(junctions)[["case_control"]] == "case"]
    assays(junctions)[["score"]] <- case_score

    print(stringr::str_c(Sys.time(), " - done!"))

    return(junctions)
}

#' Calculate z-score for x from the distribution y
#'
#' \code{.zscore} calculates a z-score for each junction count from a patient
#' sample (\code{x}), indicating it's deviation from the distribution of a
#' controls counts (\code{y}) of the same junction.
#'
#' @param x numeric vector containing case counts for 1 junction.
#' @param y numeric vector containing control counts for 1 junction.
#' @param sd_const a numeric scalar to be added to all control standard
#'   deviations. This is to prevent infinate/NaN values occuring when the
#'   standard deviation of the control counts is 0.
#'
#' @return numeric vector of length equal to \code{x} containing the z-score for
#'   the case junctions.
#'
#' @keywords internal
#' @noRd
.zscore <- function(x, y, sd_const = 0.01) {
    x_score <- (x - mean(y)) / (sd(y) + sd_const)

    return(x_score)
}
