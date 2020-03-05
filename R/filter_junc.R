#' Obtain the junctions to filter based on count
#'
#' \code{get_junc_count_filter} takes as input raw and normalised counts then
#' applies a inputted filter to obtain a lgl vector the same length as the
#' junctions determining which to keep.
#'
#' @param raw_counts df. Detailing the raw counts for each junction. Rows as
#'   junctions, cols as samples, values as raw counts.
#' @param norm_counts df. Detailing the normalised counts for each junction.
#'   Rows as junctions, cols as samples, values as raw counts.
#' @param raw_thresh dbl scalar. How many raw counts needed to keep a filter in?
#' @param norm_thresh dbl scalar. How many normalised counts needed to keep a
#'   filter in?
#' @param n_samp int scalar. Number of samples per junc that must have counts ab
#'   \code{raw_thresh} AND \code{norm_thresh}
#'
#' @return lgl vector. Each element corresponds to a junc. TRUE meaning the
#'   junction has passed the count filter.
#' @export
#'
#' @examples
get_junc_count_filter <- function(raw_counts, norm_counts, raw_thresh = 3, norm_thresh = 0.05, n_samp = 1){

  .get_count_ab <- function(x, count_thresh){

    return(sum(x >= count_thresh))

  }

  raw_count_ab_thresh <-
    raw_counts %>%
    apply(MARGIN = 1, FUN = .get_count_ab, count_thresh = raw_thresh) %>%
    unlist()

  norm_count_ab_thresh <-
    norm_counts %>%
    apply(MARGIN = 1, FUN = .get_count_ab, count_thresh = norm_thresh) %>%
    unlist()

  junc_count_filter <- (raw_count_ab_thresh >= n_samp) & (norm_count_ab_thresh >= n_samp)

  return(junc_count_filter)

}


#' Filter junctions using count and acceptor/donor annotation
#'
#' \code{filter_junc} takes as input raw and normalised counts, then applies
#' both a count and acceptor/donor filter.
#'
#' @inheritParams get_junc_count_filter
#' @param junc_count_filter lgl vector. Returned from
#'   \code{get_junc_count_filter}.
#' @param junc_metadata df. Detailing the metadata of the columns, must include
#'   "annot_ref" column with whether junc has an acceptor or donor site
#'   overlapping a known exon.
#' @param annot_ref_filter int scalar. Which types of annotation to exclude.
#'   Default is "none" or "ambig"
#'
#' @return list. Raw, normalised counts and metadata data filtered.
#'
#' @export
#'
#' @examples
filter_junc <- function(raw_counts, norm_counts, junc_count_filter, junc_metadata, annot_ref_filter = c("none", "ambig_gene"), size = NULL){

  if(!any("annot_ref" %in% names(GenomicRanges::mcols(junc_metadata)))) stop("No annotation for acceptor/donor found in metadata...")

  junc_filter <- (!(junc_metadata$annot_ref %in% annot_ref_filter) &
                    junc_count_filter)

  if("ambig_gene" %in% annot_ref_filter){

    # remove the juncs that are annotated by connect to 2 genes according to reference annotation
    junc_filter <- (junc_filter & (base::lengths(junc_metadata$gene_id_junc) < 2))

  }

  if(!is.null(size)){

    junc_filter <- (junc_filter & (GenomicRanges::width(junc_metadata) >= size[1]) & (GenomicRanges::width(junc_metadata) <= size[2]))

  }

  # filter for juncs that have an A/D site annotated and above count threshold
  raw_counts_filtered <- raw_counts %>% filter(junc_filter)
  norm_counts_filtered <- norm_counts %>% filter(junc_filter)

  # filter metadata to match junctions
  junc_metadata_filtered <- junc_metadata[junc_filter]

  # check the indexes still match
  stopifnot(identical(junc_metadata_filtered$index, raw_counts_filtered$index))

  junc_raw_norm_metadata_filtered <-
    list(raw = raw_counts_filtered,
         norm = norm_counts_filtered,
         metadata = junc_metadata_filtered)

  return(junc_raw_norm_metadata_filtered)

}
