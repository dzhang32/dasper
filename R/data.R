#' Random set of example junctions from two samples
#'
#' A dataset containing the example junction data for 2 samples outputted from
#' \code{\link{junction_load}}.
#'
#' @format RangedSummarizedExperiment from the
#'   \code{\link{SummarizedExperiment}} detailing the counts, co-ordinates of
#'   junctions for 2 example samples:
#' \describe{
#'   \item{assays}{matrix with counts for ~20,000 junctions and 2 samples}
#'   \item{colData}{example sample metadata}
#'   \item{rowRanges}{GRanges object describing the co-ordinates and strand of each junction}
#'   }
#' @source generated using data-raw/example_junctions.R
"example_junctions"
