# Load data ---------------------------------------------------------------

# use Genomic state to load txdb (GENCODE v31)
ref <- GenomicState::GenomicStateHub(version = "31", genome = "hg38", filetype = "TxDb")[[1]]

junctions_processed_path <- file.path(tempdir(), "junctions_processed.rda")

if (file.exists(junctions_processed_path)) {
    load(junctions_processed_path)
} else {
    junctions_processed <-
        junction_process(
            junctions_example,
            count_thresh = NULL,
            n_samp = NULL,
            ref,
            sd_const = 0.02
        )

    save(junctions_processed,
        file = junctions_processed_path
    )
}

# filter junctions to save time
junctions <- junctions_processed

# simulate random coverage scores to save time
coverage_scores <- matrix(
    data = sample(-100:100, dim(junctions)[1] * dim(junctions)[2], replace = TRUE),
    nrow = dim(junctions)[1],
    ncol = dim(junctions)[2]
)

colnames(coverage_scores) <- dimnames(junctions)[[2]]
assays(junctions)[["coverage_score"]] <- coverage_scores

# Sashimi plot functions --------------------------------------------------

# set one of the samples as control for testing plotting of controls
SummarizedExperiment::colData(junctions_processed)[["case_control"]] <- "control"
junctions_processed <-
    junctions_processed %>%
    junction_filter(count_thresh = c("raw" = 3))

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

# pick the gene with the most number of junctions
gene_id_to_plot <-
    SummarizedExperiment::rowData(junctions_processed)[["gene_id_junction"]] %>%
    unlist() %>%
    table() %>%
    sort() %>%
    tail(1) %>%
    names()

gene_tx_list <- .gene_tx_type_get(gene_id_to_plot)

exons_to_plot <- .exons_to_plot_get(
    ref = ref,
    gene_tx_list = gene_tx_list,
    region = NULL
)

# take half the gene as the region of interest
region <- GRanges(
    paste0(
        seqnames(exons_to_plot) %>% as.character() %>% unique(),
        ":",
        start(exons_to_plot) %>% min(),
        "-",
        mean(c(min(start(exons_to_plot)), max(end(exons_to_plot))))
    )
)

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
# pick the gene with the most number of junctions
tx_name_to_plot <-
    SummarizedExperiment::rowData(junctions_processed)[["tx_name_start"]] %>%
    unlist() %>%
    table() %>%
    sort() %>%
    tail(1) %>%
    names()

junctions_to_plot <- .junctions_to_plot_get(junctions_processed, gene_tx_list, region)

test_that("junctions_to_plot has the correct output", {
    expect_true(all(any(mcols(junctions_to_plot)[["gene_id_start"]] == gene_tx_list) |
        any(mcols(junctions_to_plot)[["gene_id_end"]] == gene_tx_list)))

    expect_equal(
        findOverlaps(junctions_to_plot, region) %>% length(),
        length(junctions_to_plot)
    )

    junctions_to_plot_2 <- .junctions_to_plot_get(junctions_processed, list(tx_name = tx_name_to_plot), region = NULL)

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
    sum_func = mean,
    digits = 3,
    assay_name = "norm"
)

junctions_to_plot_2 <- .junctions_counts_type_get(junctions_to_plot,
    case_id = list(samp_id = c("samp_1", "samp_2")),
    sum_func = NULL,
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
    junctions = junctions_processed,
    ref = ref,
    gene_tx_id = gene_id_to_plot,
    region = region,
    count_label = FALSE
)

sashimi2 <- plot_sashimi(
    junctions = junctions_processed,
    ref = ref,
    gene_tx_id = tx_name_to_plot,
    case_id = list(samp_id = c("samp_1")),
    sum_func = mean,
    count_label = TRUE,
    digits = 2
)

test_that("coords_to_plot has the correct output", {
    expect_true(is(sashimi1, "ggplot"))
    expect_true(is(sashimi2, "ggplot"))
})
