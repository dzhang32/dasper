#' Summarise outlier info per junction to cluster-level
#'
#' \code{outlier_aggregate} will aggregate the outlier scores into a
#' cluster-level. It will then rank each cluster based on this aggregated score
#' and annotate each cluster with it's most likely associated gene and
#' transcript.
#'
#' @inheritParams junction_annot
#'
#' @return a data.frame with one row per cluster detailing each cluster's
#'   associated junctions, outlier scores, genes and transcripts of each
#'   patient.
#'
#' @keywords internal
#' @noRd
outlier_aggregate <- function(junctions, samp_id_col = "samp_id") {

    ##### Check user input #####

    if (!all(c("direction", "outlier_score") %in% names(assays(junctions)))) {
        stop("junctions must contain both 'direction' and 'outlier_score' assays")
    }

    ##### Get outlier scores per sample #####

    outlier_scores_samp <- vector(mode = "list", length = dim(junctions)[2])

    for (i in seq_len(dim(junctions)[2])) {
        outlier_scores_samp[[i]] <-
            dplyr::tibble(
                samp_id = colData(junctions)[[samp_id_col]][i],
                direction = assays(junctions)[["direction"]][, i],
                outlier_score = assays(junctions)[["outlier_score"]][, i]
            )
    }

    outlier_scores_samp <- do.call(dplyr::bind_rows, outlier_scores_samp)

    ##### Aggregate scores to a cluster-level #####

    outlier_cluster



    names(clusters) <- 1:nrow(outlier_score)
    clusters <- clusters %>% unlist()

    # there's 1 seed junction per cluster
    # obtain the most disregulated UJ and DJ
    # per cluster
    cluster_outlier_score <-
        dplyr::tibble(
            cluster = names(clusters),
            index = unname(clusters)
        ) %>%
        dplyr::left_join(outlier_score) %>%
        dplyr::group_by(cluster, up_down) %>%
        dplyr::filter(
            !duplicated(outlier_score), # remove duplicated most-dysregulated scores
            outlier_score == min(outlier_score)
        ) %>%
        dplyr::ungroup()

    # the n of junctions per cluster should be 1 or 2
    # keep only clusters that have at least 1 UJ and DJ
    cluster_outlier_score <- cluster_outlier_score %>%
        dplyr::group_by(cluster) %>%
        dplyr::filter(dplyr::n() == 2) %>%
        ungroup()

    # obtain mean outlier score between most dysregulated UJ/DJ
    # str_c to obtain all the
    cluster_outlier_score <- cluster_outlier_score %>%
        dplyr::arrange(cluster, index) %>% # makes sure indexes are in same order for different clusters
        group_by(cluster) %>%
        dplyr::summarise(
            indexes = stringr::str_c(index, collapse = ";"),
            mean_outlier_score = mean(outlier_score)
        )

    # remove clusters that have identical most dysregulated UJ/DJ
    cluster_outlier_score <-
        cluster_outlier_score %>%
        dplyr::filter(!duplicated(indexes)) %>%
        dplyr::mutate(rank = rank(mean_outlier_score))

    samp_of_interest <-
        which(SummarizedExperiment::colData(junctions)[["samp_id_tidy"]] ==
            unique(outlier_score[["samp_id_tidy"]]))

    samp_of_interest <- SummarizedExperiment::colData(junctions)[samp_of_interest, ]

    which_patho <- all(gene_id_junction[cluster_outlier_score[["cluster"]] %>% as.integer()] ==
        samp_of_interest[["gene_id"]])

    cluster_outlier_score_patho <- cluster_outlier_score[which_patho, ] %>%
        mutate(samp_id_tidy = samp_of_interest[["samp_id_tidy"]])

    return(cluster_outlier_score_patho)
}

.outlier_cluster <- function(junctions, outlier_scores_samp) {

    # there's 1 seed junction per cluster
    # obtain the most disregulated UJ and DJ
    # per cluster
    cluster_outlier_score <-
        dplyr::tibble(
            cluster = names(clusters),
            index = unname(clusters)
        ) %>%
        dplyr::left_join(outlier_score) %>%
        dplyr::group_by(cluster, up_down) %>%
        dplyr::filter(
            !duplicated(outlier_score), # remove duplicated most-dysregulated scores
            outlier_score == min(outlier_score)
        ) %>%
        dplyr::ungroup()

    # the n of junctions per cluster should be 1 or 2
    # keep only clusters that have at least 1 UJ and DJ
    cluster_outlier_score <- cluster_outlier_score %>%
        dplyr::group_by(cluster) %>%
        dplyr::filter(dplyr::n() == 2) %>%
        ungroup()

    # obtain mean outlier score between most dysregulated UJ/DJ
    # str_c to obtain all the
    cluster_outlier_score <- cluster_outlier_score %>%
        dplyr::arrange(cluster, index) %>% # makes sure indexes are in same order for different clusters
        group_by(cluster) %>%
        dplyr::summarise(
            indexes = stringr::str_c(index, collapse = ";"),
            mean_outlier_score = mean(outlier_score)
        )

    # remove clusters that have identical most dysregulated UJ/DJ
    cluster_outlier_score <-
        cluster_outlier_score %>%
        dplyr::filter(!duplicated(indexes)) %>%
        dplyr::mutate(rank = rank(mean_outlier_score))

    samp_of_interest <-
        which(SummarizedExperiment::colData(junctions)[["samp_id_tidy"]] ==
            unique(outlier_score[["samp_id_tidy"]]))

    samp_of_interest <- SummarizedExperiment::colData(junctions)[samp_of_interest, ]

    which_patho <- all(gene_id_junction[cluster_outlier_score[["cluster"]] %>% as.integer()] ==
        samp_of_interest[["gene_id"]])

    cluster_outlier_score_patho <- cluster_outlier_score[which_patho, ] %>%
        mutate(samp_id_tidy = samp_of_interest[["samp_id_tidy"]])
}
