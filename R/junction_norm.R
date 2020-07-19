#' Normalise junction counts
#'
#' \code{junction_norm} normalises the junction counts for each sample. First,
#' each junction is clustered by finding all other junctions that share with it
#' an acceptor or donor position. Then, a proportion-spliced-in (PSI) is
#' calculated for each junction by dividing the raw junction count by the total
#' number of counts in it's associated cluster.
#'
#' @inheritParams junction_annot
#'
#' @return junctions as a with additional an
#'   \code{\link[SummarizedExperiment]{assay}} containing the normalised
#'   counts.
#'
#' @examples
#'
#' junctions <- junction_norm(junctions_example)
#' junctions
#' @export
junction_norm <- function(junctions) {

    ##### Cluster junctions by matching start and end #####

    print(stringr::str_c(Sys.time(), (" - Clustering junctions...")))

    junctions <- .junction_cluster(junctions)

    ##### Normalise counts #####

    print(stringr::str_c(Sys.time(), (" - Normalising junction counts...")))

    junctions <- .junction_norm_count(junctions)

    print(stringr::str_c(Sys.time(), (" - done!")))

    return(junctions)
}

#' Obtain junction clusters
#'
#' \code{.get_junc_cluster} defines a cluster for each junction. Each junction
#' is used as a seed to find other junctions that share it's acceptor or donor
#' position. Each cluster is comprised by it's seed junction along with any
#' junctions matching with a acceptor/donor position.
#'
#' @inheritParams junction_annot
#'
#' @return junctions with cluster definitions.
#'
#' @keywords internal
#' @noRd
.junction_cluster <- function(junctions) {

    # find overlaps between the junction starts and ends
    junctions_start_end <- .get_start_end(junctions)

    # only find hits between start/start or end/end
    # the other overlap (start/end) does not capture splicing events that
    # necessarily tag the same intron
    start_hits <- findOverlaps(
        query = junctions_start_end[["start"]],
        subject = junctions_start_end[["start"]],
        type = "equal",
        ignore.strand = F
    )

    end_hits <- findOverlaps(
        query = junctions_start_end[["end"]],
        subject = junctions_start_end[["end"]],
        type = "equal",
        ignore.strand = F
    )

    # query hits correspond to the index of each seed junction
    # subject hits correspond to every junction that matches each seed junctions start/end
    # this regroups subject junctions by corresponding query seed junction group
    # outputs a list containing an element for each group/cluster/seed junction
    # each element contains all junction indexes in that cluster
    clusters <- .regroup(
        x = c(subjectHits(start_hits), subjectHits(end_hits)),
        groups = c(queryHits(start_hits), queryHits(end_hits)),
        all_groups = 1:length(junctions)
    ) %>%
        IRanges::IntegerList() %>%
        unique() %>% # unique so that each junction only appears once per cluster
        sort()

    # add clusters and indexes to coords
    mcols(junctions)[["index"]] <- 1:length(junctions)
    mcols(junctions)[["clusters"]] <- clusters

    return(junctions)
}

#' Normalise junction counts
#'
#' \code{..junction_norm_count} normalises the junction counts for each sample
#' by dividing the counts for each junction by the total number of counts in
#' it's associated cluster.
#'
#' @inheritParams junction_annot
#'
#' @return junctions with an additional
#'   \code{\link[SummarizedExperiment]{assay}} containing normalised counts.
#'
#' @keywords internal
#' @noRd
.junction_norm_count <- function(junctions) {
    n_junctions <- nrow(SummarizedExperiment::assays(junctions)[["raw"]])
    n_samps <- ncol(SummarizedExperiment::assays(junctions)[["raw"]])

    norm_counts <- matrix(
        nrow = n_junctions,
        ncol = n_samps
    )

    # unlisted values are the junction indexes in each cluster
    # unlisted names are the seed junction indexes/groups of each cluster
    clusters_unlisted <- mcols(junctions)[["clusters"]] %>%
        unlist()

    # for every sample
    for (i in 1:n_samps) {

        # 1. index the raw counts by the contents of every cluster
        # 2. then split() this into a list with each element as a vector
        # containing the raw reads for every junction in that cluster
        # 3. sum over each element obtain the total counts per cluster
        cluster_sum_counts <-
            .regroup(
                x = SummarizedExperiment::assays(junctions)[["raw"]][, i][clusters_unlisted],
                group = names(clusters_unlisted),
                all_groups = 1:length(junctions)
            ) %>%
            lapply(FUN = sum) %>%
            unlist()

        # clusters may be returned in a different order
        # since .regroup()/split() converts names to character then sorts
        # re-order by the original order
        cluster_sum_counts <- cluster_sum_counts[as.character(mcols(junctions)[["index"]])]

        norm_counts[, i] <- SummarizedExperiment::assays(junctions)[["raw"]][, i] / cluster_sum_counts
    }

    colnames(norm_counts) <- colnames(SummarizedExperiment::assays(junctions)[["raw"]])

    # convert NA's (junctions that have 0 reads in their cluster) to 0
    norm_counts[is.na(norm_counts)] <- 0

    SummarizedExperiment::assays(junctions)[["norm"]] <- norm_counts

    return(junctions)
}
