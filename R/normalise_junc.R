#' Normalise junction counts
#'
#' \code{normalise_junc} normalises the junction counts for each sample. First,
#' each junction is clustered by finding all other junctions that share with it
#' an acceptor or donor position. Then, a proportion-spliced-in is calculated
#' for each junction by dividing the raw junction count by the total number of
#' counts in it's associated cluster.
#'
#' @inheritParams annotate_junc_ref
#'
#' @return list containing metadata in a \code{\link[GenomicRanges]{GRanges}}
#'   format and dataframes detailing raw and normalised counts.
#'
#' @export
normalise_junc <- function(junc_metadata, raw_count){

  ##### Cluster junctions by matching start and end #####

  print(stringr::str_c(Sys.time(), (" - Clustering junctions...")))

  junc_metadata <- .get_junc_cluster(junc_metadata)

  ##### Cluster junctions by matching start and end #####

  print(stringr::str_c(Sys.time(), (" - Normalising junction counts...")))

  norm_count <- .normalise_count(junc_metadata, raw_count)

  print(stringr::str_c(Sys.time(), (" - done!")))

  juncs <- list(metadata = junc_metadata,
                raw_count = raw_count,
                norm_count = norm_count)

  return(juncs)

}

#' Obtain junction clusters
#'
#' \code{.get_junc_cluster} defines a cluster for each junction. Each junction
#' is used as a seed to find other junctions that share it's acceptor or donor
#' position. Each cluster is comprised by it's seed junction along with any
#' junctions matching with a acceptor/donor position.
#'
#' @inheritParams annotate_junc_ref
#'
#' @return junction metadata with cluster definitions.
.get_junc_cluster <- function(junc_metadata){

  # find overlaps between the junction starts and ends
  junc_start_end <- .get_gr_for_start_end(junc_metadata)

  # only find hits between start and start of junctions or end and end
  # the other overlap (start to end) does not capture splicing events that
  # necessarily tag the same intron
  start_hits <- findOverlaps(query = junc_start_end[["start"]],
                             subject = junc_start_end[["start"]],
                             type = "equal")

  end_hits <- findOverlaps(query = junc_start_end[["end"]],
                           subject = junc_start_end[["end"]],
                           type = "equal")

  # each cluster is defined by a seed junction
  # query hits correspond to the index of the seed junction
  # subject hits correspond to every other junction that match the seed junctions start/end
  # then split() creates a list where each element is a vector
  # containing all junction indexes in that cluster
  clusters <-
    split(x = c(subjectHits(start_hits), subjectHits(end_hits)),
          f = c(queryHits(start_hits), queryHits(end_hits))) %>%
    IRanges::IntegerList() %>% # convert an IntegerList for consistency and easier manipulation
    unique() # unique so that each junction only appears once per cluster

  # add clusters and indexes to coords
  mcols(junc_metadata)[["index"]] <- 1:length(junc_metadata)
  mcols(junc_metadata)[["clusters"]] <- clusters

  return(junc_metadata)

}

#' Normalise junction counts
#'
#' \code{.normalise_count} normalises the junction counts for each sample by
#' dividing the counts for each junction by the total number of counts of it's
#' associated cluster.
#'
#' @inheritParams annotate_junc_ref
#'
#' @return dataframe matching the \code{raw_count} with columns as samples, rows
#'   as junctions with elements replaced with the normalised counts.
.normalise_count <- function(junc_metadata, raw_count){

  norm_count <- matrix(ncol = ncol(raw_count),
                       nrow = nrow(raw_count))

  # unlisted values are the junction indexes in each cluster
  # unlisted names are the seed junction indexes per cluster
  clusters_unlisted <- unlist(mcols(junc_metadata)[["clusters"]])

  # for every sample
  for(i in 1:ncol(raw_count)){

    # 1. index the raw counts by the contents of every cluster
    # 2. then split() this into a list with each element as a vector
    # containing the raw reads for every junction in that cluster
    # 3. sum over each element obtain the total counts per cluster
    cluster_sum_counts <-
      raw_count[[i]][clusters_unlisted] %>%
      S4Vectors::split(f = names(clusters_unlisted)) %>%
      lapply(FUN = sum) %>%
      unlist()

    # clusters are returned in a different order
    # since split() converts names to character then sorts
    # re-order by the original order
    cluster_sum_counts <- cluster_sum_counts[as.character(mcols(junc_metadata)[["index"]])]

    norm_count[,i] <- raw_count[[i]]/cluster_sum_counts

  }

  colnames(norm_count) <- colnames(raw_count)

  norm_count <- norm_count %>%
    dplyr::as_tibble()

  # convert NA's (junctions that have 0 reads in their cluster) to 0
  norm_count[is.na(norm_count)] <- 0

  return(norm_count)

}
