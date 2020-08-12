context("Testing outlier processing")

##### reticulate set-up #####

if (.Platform$OS.type != "windows") {

    # force reticulate to use the python3 install
    # if windows skip this step
    reticulate::use_python(Sys.which("python3"), required = TRUE)
}

##### Set up random scores data #####

# needs annotation and normalisation for
# aggregation of outlier scores
suppressWarnings(expr = {
    junctions <- junctions_example[, colData(junctions_example)[["case_control"]] == "case"] %>%
        junction_annot(ref = "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz") %>%
        junction_norm()
})

# add random scoress and direction to save time

direction <- matrix(
    data = sample(c(1, -1), dim(junctions)[1] * dim(junctions)[2], replace = TRUE),
    nrow = dim(junctions)[1],
    ncol = dim(junctions)[2]
)

scores <- matrix(
    data = sample(-100:100, dim(junctions)[1] * dim(junctions)[2], replace = TRUE),
    nrow = dim(junctions)[1],
    ncol = dim(junctions)[2]
)

coverage_scores <- matrix(
    data = sample(-100:100, dim(junctions)[1] * dim(junctions)[2], replace = TRUE),
    nrow = dim(junctions)[1],
    ncol = dim(junctions)[2]
)

colnames(direction) <- dimnames(junctions)[[2]]
colnames(scores) <- dimnames(junctions)[[2]]
colnames(coverage_scores) <- dimnames(junctions)[[2]]

assays(junctions)[["direction"]] <- direction
assays(junctions)[["score"]] <- scores
assays(junctions)[["coverage_score"]] <- coverage_scores

# Detecting outliers using an isolation forest ----------------------------

##### outlier_detect #####

junctions_w_outlier_scores <- outlier_detect(junctions,
    feature_names = c("score", "coverage_score"),
    random_state = 32L
)

up_indexes <- which(assays(junctions)[["direction"]][, 1] == 1)
outlier_up <- .outlier_score(
    features =
        data.frame(
            scores = assays(junctions)[["score"]][, 1][up_indexes],
            coverage_scores = assays(junctions)[["coverage_score"]][, 1][up_indexes]
        ),
    random_state = 32L
)

down_indexes <- which(assays(junctions)[["direction"]][, 1] == -1)
outlier_down <- .outlier_score(
    features =
        data.frame(
            scores = assays(junctions)[["score"]][, 1][down_indexes],
            coverage_scores = assays(junctions)[["coverage_score"]][, 1][down_indexes]
        ),
    random_state = 32L
)

test_that("outlier_detect has the correct output", {

    # if junctions have been reordered either one of the direction
    # or scoress of up/down would not match
    expect_identical(
        assays(junctions_w_outlier_scores)[["direction"]],
        assays(junctions)[["direction"]]
    )

    expect_identical(
        outlier_up %>% as.numeric(),
        assays(junctions_w_outlier_scores)[["outlier_score"]][, 1][up_indexes]
    )

    expect_identical(
        outlier_down %>% as.numeric(),
        assays(junctions_w_outlier_scores)[["outlier_score"]][, 1][down_indexes]
    )
})

test_that("outlier_detect catches user-input errors", {
    expect_error(
        outlier_detect(junctions, feature_names = c("not_an_assay")),
        "Assays does not contain the following: "
    )

    assays(junctions)[["direction"]] <- NULL

    expect_error(
        outlier_detect(junctions),
        "junctions must contain a 'direction' assay"
    )
})

# Aggregating outlier scores to cluster-level -----------------------------

##### .outlier_wrangle #####

outlier_scores_samp <- .outlier_wrangle(junctions_w_outlier_scores, samp_id_col = "samp_id")

test_that(".outlier_wrangle has the correct output", {
    expect_true(is(outlier_scores_samp, "list"))
    expect_true(all(lapply(outlier_scores_samp, nrow) == length(junctions_w_outlier_scores)))

    expect_identical(
        outlier_scores_samp[[1]][["direction"]],
        assays(junctions_w_outlier_scores)[["direction"]][, 1]
    )
    expect_identical(
        outlier_scores_samp[[2]][["direction"]],
        assays(junctions_w_outlier_scores)[["direction"]][, 2]
    )

    expect_identical(
        outlier_scores_samp[[1]][["outlier_score"]],
        assays(junctions_w_outlier_scores)[["outlier_score"]][, 1]
    )
    expect_identical(
        outlier_scores_samp[[2]][["outlier_score"]],
        assays(junctions_w_outlier_scores)[["outlier_score"]][, 2]
    )
})

##### .outlier_cluster #####

outlier_scores_samp <- .outlier_cluster(junctions_w_outlier_scores, outlier_scores_samp)

outlier_cluster_check <- function(junctions_w_outlier_scores, outlier_scores_samp) {
    clusters <- unlist(SummarizedExperiment::rowData(junctions_w_outlier_scores)[["clusters"]])

    junctions_by_cluster <- junctions_w_outlier_scores[unname(clusters)]

    check <- TRUE

    for (i in seq_along(outlier_scores_samp)) {
        outlier_clusters <-
            lapply(assays(junctions_by_cluster), FUN = function(x) x[, i]) %>%
            as.data.frame() %>%
            dplyr::mutate(cluster_index = names(clusters))

        outlier_clusters <- outlier_clusters %>%
            dplyr::group_by(cluster_index, direction) %>%
            dplyr::filter(
                !duplicated(outlier_score),
                outlier_score == min(outlier_score)
            ) %>%
            dplyr::group_by(cluster_index) %>%
            dplyr::filter(dplyr::n() != 1)

        outlier_clusters <- outlier_clusters %>%
            dplyr::group_by(cluster_index) %>%
            dplyr::summarise(mean_outlier_score = mean(outlier_score))

        check <- all(check, all(outlier_scores_samp[[i]][["cluster_index"]] %in%
            outlier_clusters[["cluster_index"]]))

        check <- all(check, identical(
            sort(unique(outlier_scores_samp[[i]][["mean_outlier_score"]])),
            sort(unique(outlier_clusters[["mean_outlier_score"]]))
        ))
    }

    return(check)
}

test_that(".outlier_cluster has the correct output", {
    expect_true(is(outlier_scores_samp, "list"))
    expect_true(all(unlist(lapply(outlier_scores_samp[[1]], nrow)) == 2))
    expect_true(all(unlist(lapply(outlier_scores_samp[[2]], nrow)) == 2))

    expect_true(outlier_cluster_check(
        junctions_w_outlier_scores,
        outlier_scores_samp
    ))
})

##### .outlier_cluster_tidy #####

outlier_scores_tidy <- .outlier_cluster_tidy(outlier_scores_samp)

test_that(".outlier_cluster_tidy has the correct output", {
    expect_true(is(outlier_scores_tidy, "DataFrame"))
    expect_true(is(outlier_scores_tidy[["gene_id_cluster"]], "CharacterList"))
    expect_true(is(outlier_scores_tidy[["junctions"]], "list"))
    expect_true(all(colData(junctions_w_outlier_scores)[["samp_id"]] %in%
        unique(outlier_scores_tidy[["samp_id"]])))
})

##### outlier_aggregate #####

test_that("outlier_aggregate catches user input errors", {
    expect_error(
        outlier_aggregate(junctions_example),
        "junctions rowData must contain 'clusters'. Have you run junction_norm?"
    )

    expect_error(
        outlier_aggregate(junctions),
        "junctions must contain both 'direction' and 'outlier_score' assays"
    )

    SummarizedExperiment::rowData(junctions_w_outlier_scores)[["clusters"]][[1]] <-
        length(junctions_w_outlier_scores) + 1

    expect_error(
        outlier_aggregate(junctions_w_outlier_scores),
        "Not all cluster indexes match junctions. Have you filtered junctions after running junction_norm?"
    )
})

# Process outliers --------------------------------------------------------

outlier_scores_tidy_2 <- outlier_process(junctions, random_state = 32L)

test_that(".outlier_cluster_tidy has the correct output", {
    expect_equivalent(outlier_scores_tidy, outlier_scores_tidy_2)
})
