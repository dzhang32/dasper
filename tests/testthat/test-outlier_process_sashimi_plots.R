context("Testing outlier processing and sashimi plot functions")

##### Set up random scores data #####

# use Genomic state to load txdb (GENCODE v31)
ref <- GenomicState::GenomicStateHub(version = "31", genome = "hg38", filetype = "TxDb")[[1]]

# needs annotation and normalisation for
# aggregation of outlier scores
# and to test sashimi plotting
junctions_annoted_normed <- junctions_example %>%
    junction_filter() %>%
    junction_annot(ref) %>%
    junction_filter(types = c("unannotated", "ambig_gene")) %>%
    junction_norm()

# take first 100 junctions to save time
# remove cases to simulate junction_score
# for testing outlier scores
junctions <- junctions_annoted_normed[1:100, colData(junctions_annoted_normed)[["case_control"]] == "case"]

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

# Generate key of cluster indexes vs junction indexes
clusters <- SummarizedExperiment::rowData(junctions)[["clusters"]]
names(clusters) <- seq_len(dim(junctions)[1])
clusters <- unlist(clusters)
clusters <- dplyr::tibble(
    cluster_index = names(clusters),
    junction_index = unname(clusters)
)

outlier_scores_samp <- BiocParallel::bplapply(outlier_scores_samp,
    FUN = .outlier_cluster,
    BPPARAM = BiocParallel::SerialParam(),
    clusters = clusters
)

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

##### .outlier_cluster_annot & .outlier_cluster_tidy #####

outlier_scores_samp <- BiocParallel::bplapply(outlier_scores_samp,
    FUN = .outlier_cluster_annot,
    BPPARAM = BiocParallel::SerialParam()
)

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

outlier_scores_tidy_2 <- outlier_process(junctions,
    random_state = 32L
)

test_that(".outlier_cluster_tidy has the correct output", {
    expect_equivalent(outlier_scores_tidy, outlier_scores_tidy_2)
})

# Sashimi plot functions --------------------------------------------------

SummarizedExperiment::rowData(junctions_annoted_normed)[["tx_name_start"]] %>%
    unlist() %>%
    table() %>%
    sort() %>%
    tail()

##### .gene_tx_type_get #####

test_that(".gene_tx_type_get has the correct output", {
    expect_equal(
        .gene_tx_type_get("ENSG_example"),
        list(gene_id = "ENSG_example")
    )

    expect_equal(
        .gene_tx_type_get("ENST_example"),
        list(tx_name = "ENST_example")
    )
})

test_that(".gene_tx_type_get catches user-input errors", {
    expect_error(
        .gene_tx_type_get(NULL),
        "gene_tx_id must be set and be of length 1"
    )

    expect_error(
        .gene_tx_type_get(c("ENSG_example", "ENSG_example")),
        "gene_tx_id must be set and be of length 1"
    )

    # for now, dasper does not accept gene symbols
    expect_error(
        .gene_tx_type_get(c("SNCA")),
        "gene_tx_id does not include an ENST or ENSG prefix"
    )
})

##### .exons_to_plot_get #####

gene_id_to_plot <- "ENSG00000241973.10"
gene_tx_list <- .gene_tx_type_get(gene_id_to_plot)
region <- GRanges("chr22:20740000-20760000")

exons_to_plot <- .exons_to_plot_get(
    ref = ref,
    gene_tx_list = gene_tx_list,
    region = region
)

test_that(".exons_to_plot_get has the correct output", {
    expect_equal(
        .exons_to_plot_get(ref,
            .gene_tx_type_get(gene_id_to_plot),
            region = NULL
        ),
        .gene_tx_type_get(gene_id_to_plot) %>%
            GenomicFeatures::exons(ref, filter = .) %>%
            GenomicRanges::disjoin()
    )

    expect_equal(
        findOverlaps(exons_to_plot, region) %>% length(),
        length(exons_to_plot)
    )
})

test_that(".junctions_to_plot_get catches user input errors", {
    expect_error(
        .exons_to_plot_get(ref, gene_tx_list,
            region = GRanges(c("chr22:1-1", "chr22:1-1"))
        ),
        "region must be a GenomicRanges object of length 1"
    )

    expect_error(
        .exons_to_plot_get(ref, gene_tx_list,
            region = data.frame()
        ),
        "region must be a GenomicRanges object of length 1"
    )

    expect_error(
        .exons_to_plot_get(ref, list(gene_id = "ENSG_not_a_real_gene"),
            region = NULL
        ),
        "No exons found to plot"
    )

    expect_error(
        .exons_to_plot_get(ref, gene_tx_list,
            region = GRanges("chr22:1-1")
        ),
        "No exons found to plot"
    )
})

##### .junctions_to_plot_get #####

# select a transcript of interest
# mcols(junctions)[["tx_name_start"]] %>% unlist() %>% table() %>% sort(decreasing = TRUE) %>% head()
tx_name_to_plot <- "ENST00000255882.10"

junctions_to_plot <- .junctions_to_plot_get(junctions_annoted_normed, gene_tx_list, region)

test_that("junctions_to_plot has the correct output", {
    expect_true(all(any(mcols(junctions_to_plot)[["gene_id_start"]] == gene_tx_list) |
        any(mcols(junctions_to_plot)[["gene_id_end"]] == gene_tx_list)))

    expect_equal(
        findOverlaps(junctions_to_plot, region) %>% length(),
        length(junctions_to_plot)
    )

    junctions_to_plot_2 <- .junctions_to_plot_get(junctions_annoted_normed, list(tx_name = tx_name_to_plot), region = NULL)

    expect_true(all(any(mcols(junctions_to_plot_2)[["tx_name_start"]] == tx_name_to_plot) |
        any(mcols(junctions_to_plot_2)[["tx_name_end"]] == tx_name_to_plot)))
})

test_that(".junctions_to_plot_get catches user input errors", {
    expect_error(
        .junctions_to_plot_get(junctions,
            list(gene_id = "ENSG_not_a_real_gene"),
            region = NULL
        ),
        "No junctions found to plot"
    )
})

##### .coords_to_plot_get #####

gene_tx_to_plot <- GenomicFeatures::genes(ref, filter = gene_tx_list)
coords_to_plot <- .coords_to_plot_get(gene_tx_to_plot, exons_to_plot, junctions_to_plot, ext_factor = 20)

test_that("coords_to_plot has the correct output", {
    expect_equal(
        coords_to_plot[["chr"]],
        seqnames(exons_to_plot) %>%
            as.character() %>% unique()
    )

    expect_equal(
        coords_to_plot[["strand"]],
        strand(exons_to_plot) %>%
            as.character() %>% unique()
    )

    min_start <- min(c(start(exons_to_plot), start(junctions_to_plot)))
    max_end <- max(c(end(exons_to_plot), end(junctions_to_plot)))

    expect_equal(
        coords_to_plot[["x_min"]],
        min_start - (max_end - min_start) / 20
    )

    expect_equal(
        coords_to_plot[["x_max"]],
        max_end + (max_end - min_start) / 20
    )
})

##### .plot_gene_track #####

gene_track <- .plot_gene_track(coords_to_plot, exons_to_plot)

test_that("coords_to_plot has the correct output", {
    expect_true(is(gene_track, "ggplot"))
})

##### .junctions_counts_type_get #####

junctions_to_plot_1 <- .junctions_counts_type_get(junctions_to_plot,
    case_id = list(samp_id = "samp_1"),
    control_agg_func = mean,
    digits = 3,
    assay_name = "norm"
)

junctions_to_plot_2 <- .junctions_counts_type_get(junctions_to_plot,
    case_id = list(samp_id = c("samp_1", "samp_2")),
    control_agg_func = NULL,
    digits = 0,
    assay_name = "raw"
)

test_that("coords_to_plot has the correct output", {
    expect_true(is(junctions_to_plot_1, "data.frame"))
    expect_true(is(junctions_to_plot_2, "data.frame"))
    expect_true(length(junctions_to_plot_1) > 1)
    expect_true(length(junctions_to_plot_2) > 1)

    expect_true(all(c("samp_1", "control") %in% colnames(junctions_to_plot_1)))
    expect_true(all(c("samp_1", "samp_2") %in% colnames(junctions_to_plot_2)))
    expect_false("control" %in% colnames(junctions_to_plot_2))

    check_count_val <- function(junctions_to_plot, samp_ids, range) {
        counts <- junctions_to_plot[, samp_ids] %>%
            unlist()

        all(counts >= range[1] & counts <= range[2])
    }

    expect_true(check_count_val(
        junctions_to_plot_1,
        c("samp_1", "control"),
        c(0, 1)
    ))

    expect_false(check_count_val(
        junctions_to_plot_2,
        c("samp_1", "samp_2"),
        c(0, 1)
    ))
})

##### plot_sashimi #####

sashimi1 <- plot_sashimi(
    junctions = junctions_annoted_normed,
    ref = ref,
    gene_tx_id = gene_id_to_plot,
    count_label = FALSE
)

sashimi2 <- plot_sashimi(
    junctions = junctions_annoted_normed,
    ref = ref,
    gene_tx_id = tx_name_to_plot,
    case_id = list(samp_id = c("samp_1")),
    control_agg_func = mean,
    region = region,
    count_label = TRUE,
    digits = 2
)

test_that("coords_to_plot has the correct output", {
    expect_true(is(sashimi1, "ggplot"))
    expect_true(is(sashimi2, "ggplot"))
})
