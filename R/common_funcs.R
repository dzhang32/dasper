#' Splits a \code{GRanges} object by it's start and end.
#'
#' \code{.get_gr_for_start_end} takes a \code{GRanges} object and generates 2,
#' one containing only the start co-ordinate and the other, the end.
#'
#' @param gr Any \code{GRanges} object.
#'
#' @return list of 2 grs, each with 1 range corresponding to every range in the
#'   input. One contains start positions, the other ends.
.get_gr_for_start_end <- function(gr){

  gr_start <- gr
  GenomicRanges::end(gr_start) <- GenomicRanges::start(gr_start)

  gr_end <- gr
  GenomicRanges::start(gr_end) <- GenomicRanges::end(gr_end)

  gr_start_end_list <- list(start = gr_start, end = gr_end)

  return(gr_start_end_list)

}
