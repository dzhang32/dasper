#' Filter junctions
#'
#' \code{filter_junc} filters out "noisy" junctions based on raw counts, the
#' width of junctions, annotation category of the junction returned from
#' \code{\link{annotate_junc_ref}} and whether the junction overlaps with a set
#' of (blacklist) regions.
#'
#' @inheritParams annotate_junc_ref
#'
#' @param raw_count dataframe with columns as samples, rows as junctions, cells
#'   as the raw count.
#' @param count_thresh the number of counts below which a junction will be
#'   filtered out.
#' @param n_samp the number of samples that have to express the junction above
#'   the \code{count_thresh}.
#' @param junc_width numeric vector of length 2. The first element denoting the
#'   lower limit of junction width and the second the upper limit. Junctions
#'   with widths outside this range will be filtered out.
#' @param blacklist_gr any junctions overlapping this set of regions (in a
#'   \code{\link[GenomicRanges]{GRanges}} format) will be filtered out.
#'
#' @return list containing the filtered set of junctions - metadata in a
#'   \code{\link[GenomicRanges]{GRanges}} format and raw counts as a dataframe.
#'
#' @export
filter_junc <- function(junc_metadata, raw_count,
                        count_thresh = NULL, n_samp = 1, junc_width = NULL, junc_cats = NULL, blacklist_gr = NULL){

  # keep in all junctions by default
  junc_filter <- !logical(length = length(junc_metadata))

  print(stringr::str_c(Sys.time(), " - Filtering junctions..."))

  if(!is.null(count_thresh)){

    count_filter <- .get_count_filter(raw_count, count_thresh, n_samp)
    junc_filter <- junc_filter & count_filter

  }

  if(!is.null(junc_width)){

    junc_filter <- junc_filter & (width(junc_metadata) >= junc_width[1]) & (width(junc_metadata) <= junc_width[2])

  }

  if(!is.null(junc_cats)){

    # check that user has inputted junction categories that are expected
    exp_junc_cats <- c("annotated", "novel_acceptor", "novel_donor", "novel_exon_skip", "novel_combo", "ambig_gene", "none")
    mismatch <- junc_cats[!(junc_cats %in% exp_junc_cats)]

    if(length(mismatch) != 0){

      warning(stringr::str_c("The following junction categories are not expected so ignored:\n",
                             stringr::str_c(mismatch, collapse = ", ")))

    }

    junc_filter <- junc_filter & !(mcols(junc_metadata)[["junc_cat"]] %in% junc_cats)

  }

  if(!is.null(blacklist_gr)){

    hits <- findOverlaps(junc_metadata, blacklist_gr)
    junc_filter <- junc_filter & !(1:length(junc_filter) %in% queryHits(hits))

  }

  print(stringr::str_c(Sys.time(), " - done!"))

  juncs_filtered <- list(metadata = junc_metadata[junc_filter],
                         raw_count = raw_count[junc_filter,])

  return(juncs_filtered)

}

#' Filter junctions
#'
#' \code{get_junc_count_filter} obtains a logical vector with TRUE marking the
#' junctions that are above the count filter specified.
#'
#' @inheritParams filter_junc
#'
#' @return logical vector with TRUE meaning the junction has passed the count
#'   filter.
.get_count_filter <- function(raw_count, count_thresh, n_samp){

  raw_count_ab_thresh <-
    raw_count %>%
    apply(MARGIN = 1,
          FUN = function(x, count_thresh) sum(x >= count_thresh),
          count_thresh = count_thresh) %>%
    unlist()

  raw_count_filter <- (raw_count_ab_thresh >= n_samp)

  return(raw_count_filter)

}

