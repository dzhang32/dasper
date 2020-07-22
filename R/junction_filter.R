#' Filter junctions
#'
#' \code{junction_filter} filters out "noisy" junctions based on counts, the
#' width of junctions, annotation category of the junction returned from
#' \code{\link{junction_annot}} and whether the junction overlaps with a set of
#' (blacklist) regions.
#'
#' @inheritParams junction_annot
#'
#' @param count_thresh named vector with names matching the names of the
#'   \code{\link[SummarizedExperiment]{assays}} in \code{junctions}. Values
#'   denote the number of counts below which a junction will be filtered out.
#' @param n_samp named vector with names matching the names of the
#'   \code{\link[SummarizedExperiment]{assays}} in \code{junctions}. Values
#'   denotes number of samples that have to express the junction above the
#'   \code{count_thresh} in order for that junction to not be filtered.
#' @param width_range numeric vector of length 2. The first element denoting the
#'   lower limit of junction width and the second the upper limit. Junctions
#'   with widths outside this range will be filtered out.
#' @param types any junctions matching these types, derived form
#'   \code{\link{junction_annot}} e.g. "unannotated" will be filtered out.
#' @param regions any junctions overlapping this set of regions (in a
#'   \code{\link[GenomicRanges]{GRanges}} format) will be filtered out.
#'
#' @return
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   containing filtered set of junctions.
#'
#' @seealso ENCODE blacklist regions recommended to  be included as
#'   \code{regions} can be downloaded from
#'   \url{https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz}.
#'    Further information can be found via the publication
#'   \url{https://www.nature.com/articles/s41598-019-45839-z}.
#'
#' @examples
#'
#' junctions <- junction_filter(junctions_example)
#' junctions
#' @export
junction_filter <- function(junctions,
    count_thresh = c("raw" = 5),
    n_samp = c("raw" = 1),
    width_range = NULL,
    types = NULL,
    regions = NULL) {

    # keep in all junctions by default
    junctions_filter <- !logical(length = length(junctions))

    print(stringr::str_c(Sys.time(), " - Filtering junctions..."))

    # by count
    if (!is.null(count_thresh)) {
        count_filter <- .junction_filter_count(junctions, count_thresh, n_samp)
        junctions_filter <- junctions_filter & count_filter
    }

    # by width
    if (!is.null(width_range)) {
        width_filter <- (width(junctions) >= width_range[1]) & (width(junctions) <= width_range[2])
        junctions_filter <- junctions_filter & width_filter
    }

    # by type
    if (!is.null(types)) {

        # check that user has inputted junction categories that are expected
        exp_types <- c(
            "annotated", "novel_acceptor", "novel_donor", "novel_exon_skip",
            "novel_combo", "ambig_gene", "unannotated"
        )

        mismatch <- exp_types[!(exp_types %in% exp_types)]

        if (length(mismatch) != 0) {
            warning(stringr::str_c(
                "The following junction categories are not expected so ignored:\n",
                stringr::str_c(mismatch, collapse = ", ")
            ))
        }

        type_filter <- !(mcols(junctions)[["type"]] %in% types)
        junctions_filter <- junctions_filter & type_filter
    }

    # by region
    if (!is.null(regions)) {
        hits <- findOverlaps(junctions, regions)
        region_filter <- !(seq_along(junctions_filter) %in% queryHits(hits))
        junctions_filter <- junctions_filter & region_filter
    }

    print(stringr::str_c(Sys.time(), " - done!"))

    junctions <- junctions[junctions_filter]

    return(junctions)
}

#' Filter junctions
#'
#' \code{.junction_filter_count} obtains a logical vector with TRUE marking the
#' junctions that are above the count filter specified.
#'
#' @inheritParams junction_filter
#'
#' @return logical vector with length equal to the number of \code{junctions}
#'   with TRUE meaning the junction has passed the count filter.
#'
#' @keywords internal
#' @noRd
.junction_filter_count <- function(junctions, count_thresh, n_samp) {
    count_filter_all <- !logical(dim(junctions)[1])

    for (i in seq_along(count_thresh)) {

        # how many samples for each junction have a count above count_thresh
        count_filter <-
            junctions %>%
            assays() %>%
            .[[names(count_thresh)[i]]] %>%
            apply(
                MARGIN = 1,
                FUN = function(x, count_thresh) {
                    sum(x >= count_thresh)
                },
                count_thresh = count_thresh[i]
            )

        count_filter <- count_filter >= n_samp

        count_filter_all <- count_filter_all & count_filter
    }

    return(count_filter_all)
}
