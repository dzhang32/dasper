#' Aggregate outlier scores from per junction to cluster-level
#'
#' `outlier_aggregate` will aggregate the outlier scores into a cluster-level.
#' It will then rank each cluster based on this aggregated score and annotate
#' each cluster with it's associated gene and transcript.
#'
#' @inheritParams junction_annot
#' @param samp_id_col name of the column in the
#'   [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class]  that details the sample ids.
#'
#' @return \code{DataFrame} with one row per cluster detailing each cluster's
#'   associated junctions, outlier scores, ranks and genes.
#'
#' @family outlier
#' @export
outlier_aggregate <- function(junctions, samp_id_col = "samp_id") {

    ##### Check user input #####

    if (!("clusters" %in% colnames(SummarizedExperiment::rowData(junctions)))) {
        stop("junctions rowData must contain 'clusters'. Have you run junction_norm?")
    }

    if (!all(c("direction", "outlier_score") %in% names(assays(junctions)))) {
        stop("junctions must contain both 'direction' and 'outlier_score' assays")
    }

    if (!all(unlist(SummarizedExperiment::rowData(junctions)[["clusters"]]) %in% seq_len(dim(junctions)[1]))) {
        stop("Not all cluster indexes match junctions. Have you filtered junctions after running junction_norm?")
    }

    ##### Get outlier scores per sample #####

    outlier_scores_samp <- .outlier_wrangle(junctions, samp_id_col)

    ##### Aggregate scores to a cluster-level #####

    outlier_scores_samp <- .outlier_cluster(junctions, outlier_scores_samp)

    ##### Tidy cluster-level scores #####

    outlier_scores_tidy <- .outlier_cluster_tidy(outlier_scores_samp)

    print(stringr::str_c(Sys.time(), " - done!"))

    return(outlier_scores_tidy)
}

#' Wrangles the outlier scores into a list
#'
#' `.outlier_wrangle` will extract the outlier scores from `junctions` other
#' necessary information in a format required for `.outlier_cluster`.
#'
#' @inheritParams junction_annot
#' @inheritParams outlier_aggregate
#'
#' @return a list an with one element per sample, each containing a data.frame.
#'
#' @keywords internal
#' @noRd
.outlier_wrangle <- function(junctions, samp_id_col) {
    outlier_scores_samp <- vector(mode = "list", length = dim(junctions)[2])

    for (i in seq_len(dim(junctions)[2])) {
        outlier_scores_samp[[i]] <-
            dplyr::tibble(
                samp_id = colData(junctions)[[samp_id_col]][i],
                direction = assays(junctions)[["direction"]][, i],
                outlier_score = assays(junctions)[["outlier_score"]][, i],
                gene_id_junction = SummarizedExperiment::rowData(junctions)[["gene_id_junction"]] %>% as.list()
            ) %>%
            dplyr::mutate(junction_index = dplyr::row_number())
    }

    return(outlier_scores_samp)
}

#' Aggregate the outlier scores from per junction into per cluster
#'
#' `.outlier_cluster` aggregates the scores from per junction into per
#' cluster. In order: 1. finds most dysregulated UJ and DJ 2. filters out
#' clusters without at least 1 UJ and DJ. 3. Obtain the mean outlier score
#' between the most dysregulated UJ/DJ. 4. Rank clusters by the mean outlier
#' score.
#'
#' @inheritParams outlier_scores_samp
#'
#' @return a list an with one element per sample, each containing a
#'   `data.frame`.
#'
#' @keywords internal
#' @noRd
.outlier_cluster <- function(junctions, outlier_scores_samp) {

    ##### Generate key of cluster indexes vs junction indexes #####

    clusters <- SummarizedExperiment::rowData(junctions)[["clusters"]]
    names(clusters) <- seq_len(dim(junctions)[1])
    clusters <- unlist(clusters)
    clusters <- dplyr::tibble(
        cluster_index = names(clusters),
        junction_index = unname(clusters)
    )

    for (i in seq_along(outlier_scores_samp)) {
        print(stringr::str_c(
            Sys.time(), " - Aggregating outlier scores to cluster level for sample ",
            i, "/", length(outlier_scores_samp), "..."
        ))

        outlier_clusters <- clusters %>%
            dplyr::left_join(outlier_scores_samp[[i]], by = "junction_index")

        # for each cluster, obtain most dysregulated UJs/DJs
        # remove duplicated scores so clusters are forced to contain
        # 1 or 2 junctions after this
        outlier_clusters <- outlier_clusters %>%
            dplyr::group_by(cluster_index, direction) %>%
            dplyr::filter(
                !duplicated(outlier_score),
                outlier_score == min(outlier_score)
            ) %>%
            dplyr::ungroup()

        # keep only clusters that have at least 1 UJ and DJ
        outlier_clusters <- outlier_clusters %>%
            dplyr::group_by(cluster_index) %>%
            dplyr::filter(dplyr::n() == 2) %>%
            dplyr::ungroup()

        # obtain mean outlier score between most dysregulated UJ/DJ
        outlier_clusters <- outlier_clusters %>%
            dplyr::arrange(cluster_index, junction_index) %>% # makes sure indexes are in same order for different clusters
            dplyr::group_by(cluster_index) %>%
            dplyr::mutate(mean_outlier_score = mean(outlier_score)) %>%
            dplyr::ungroup()

        # nest junction details into a list and
        # remove clusters with duplicated most dysregulated UJs and DJs
        outlier_clusters <- outlier_clusters %>%
            tidyr::nest(junctions = dplyr::one_of(c(
                "junction_index",
                "direction",
                "outlier_score",
                "gene_id_junction"
            ))) %>%
            dplyr::filter(!duplicated(junctions)) %>%
            dplyr::mutate(rank = rank(mean_outlier_score))

        outlier_scores_samp[[i]] <- outlier_clusters
    }

    return(outlier_scores_samp)
}

#' Tidy the cluster-level outlier scores
#'
#' `.outlier_cluster_tidy` tidies the the cluster-level outlier scores
#' returned by `.outlier_cluster` into a easily human-readable report.
#'
#' @inheritParams outlier_scores_samp
#'
#' @return `DataFrame` with each row corresponding to a cluster.
#'
#' @keywords internal
#' @noRd
.outlier_cluster_tidy <- function(outlier_scores_samp) {
    for (i in seq_along(outlier_scores_samp)) {
        print(stringr::str_c(
            Sys.time(), " - Tidying cluster-level outlier scores for sample ",
            i, "/", length(outlier_scores_samp), "..."
        ))

        gene_id_cluster <-
            outlier_scores_samp[[i]][["junctions"]] %>%
            lapply(FUN = function(x) {
                x[["gene_id_junction"]] %>%
                    unlist() %>%
                    unique()
            })

        outlier_scores_samp[[i]][["gene_id_cluster"]] <- gene_id_cluster
    }

    # bind all samples together
    # then convert to Dataframe to allow CharacterList column
    outlier_scores_tidy <- do.call(dplyr::bind_rows, outlier_scores_samp) %>%
        dplyr::arrange(samp_id, rank) %>%
        dplyr::select(
            samp_id,
            cluster_index,
            mean_outlier_score,
            rank,
            gene_id_cluster,
            junctions
        ) %>%
        DataFrame()

    outlier_scores_tidy[["gene_id_cluster"]] <- outlier_scores_tidy[["gene_id_cluster"]] %>%
        CharacterList()

    return(outlier_scores_tidy)
}
