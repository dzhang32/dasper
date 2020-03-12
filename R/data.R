#' Random set of junctions from two example samples
#'
#' A dataset containing the example junction data for 2 samples outputted from
#' \code{\link{merge_junc}} used for vignettes and testing.
#'
#' @format list detailing 19733 junctions with metadata:
#' \describe{
#'   \item{chr}{chromsome}
#'   \item{start}{first position of intron}
#'   \item{end}{last position of intron}
#'   \item{strand}{strand of intron}
#'   }
#'   and counts:
#' \describe{
#'   \item{eg1}{raw junction counts for sample "eg1"}
#'   \item{eg2}{raw junction counts for sample "eg1"}
#'   }
#' @source generated using data-raw/example_juncs.R
"example_juncs"
