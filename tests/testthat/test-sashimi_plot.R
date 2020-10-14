context("Test sashimi plot functions")

# use Genomic state to load txdb (GENCODE v31)
ref <- GenomicState::GenomicStateHub(version = "31", genome = "hg38", filetype = "TxDb")[[1]]

# convert seqlevels to match junctions
seqlevels(ref) <- seqlevels(ref) %>% stringr::str_replace("chr", "")

junctions <- junctions_example %>%
    junction_filter() %>%
    junction_annot(ref) %>%
    junction_filter(types = c("ambig_gene", "unannotated")) %>%
    junction_norm()

SummarizedExperiment::rowData(junctions)[["tx_name_start"]] %>%
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
region <- GRanges("22:20740000-20760000")

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
            region = GRanges(c("22:1-1", "22:1-1"))
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
            region = GRanges("22:1-1")
        ),
        "No exons found to plot"
    )
})

##### .junctions_to_plot_get #####

# select a transcript of interest
# mcols(junctions)[["tx_name_start"]] %>% unlist() %>% table() %>% sort(decreasing = TRUE) %>% head()
tx_name_to_plot <- "ENST00000255882.10"

junctions_to_plot <- .junctions_to_plot_get(junctions, gene_tx_list, region)

test_that("junctions_to_plot has the correct output", {
    expect_true(all(any(mcols(junctions_to_plot)[["gene_id_start"]] == gene_tx_list) |
        any(mcols(junctions_to_plot)[["gene_id_end"]] == gene_tx_list)))

    expect_equal(
        findOverlaps(junctions_to_plot, region) %>% length(),
        length(junctions_to_plot)
    )

    junctions_to_plot_2 <- .junctions_to_plot_get(junctions, list(tx_name = tx_name_to_plot), region = NULL)

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
    junctions = junctions,
    ref = ref,
    gene_tx_id = gene_id_to_plot,
    count_label = FALSE
)

sashimi2 <- plot_sashimi(
    junctions = junctions,
    ref = ref,
    gene_tx_id = gene_id_to_plot,
    case_id = list(samp_id = c("samp_1", "samp_2")),
    region = region,
    count_label = TRUE,
    digits = 0
)

sashimi3 <- plot_sashimi(
    junctions = junctions,
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
    expect_true(is(sashimi3, "ggplot"))
})
