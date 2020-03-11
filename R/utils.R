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
  end(gr_start) <- start(gr_start)

  gr_end <- gr
  start(gr_end) <- end(gr_end)

  gr_start_end_list <- list(start = gr_start,
                            end = gr_end)

  return(gr_start_end_list)

}

#' Merges two Characterlists into one with element-wise concatenation of the
#' vectors inside each list
#'
#' \code{.merge_lists}
#'
#' @param x CharacterList
#' @param y CharacterList
#'
#' @return
.merge_lists <- function(x, y){

  if(!identical(names(x), names(y))) stop("names of x and y lists should be identical!")

  x_y <- c(x, y) %>% unlist()

  x_y_merged <-
    x_y %>%
    unname() %>%
    split(f = names(x_y) %>%
                       factor(levels = names(x))) %>% # required to keep all levels/names
    CharacterList() %>%
    unique()

  return(x_y_merged)

}
