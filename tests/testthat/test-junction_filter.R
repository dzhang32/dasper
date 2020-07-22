context("Testing junction filtering")

##### junction_filter #####

# annotate junctions to test type filter
suppressWarnings(expr = {
    ref <- GenomicFeatures::makeTxDbFromGFF(
        "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
    )
})

junctions <- junction_annot(junctions_example, ref)

by_count <- junction_filter(junctions,
    count_thresh = c("raw" = 5),
    n_samp = c("raw" = 1),
    types = NULL
)

by_width <- junction_filter(junctions,
    count_thresh = NULL,
    n_samp = NULL,
    width_range = c(20, 1000),
    types = NULL
)

by_type <- junction_filter(junctions,
    count_thresh = NULL,
    n_samp = NULL,
    types = c("ambig_gene", "unannotated")
)

by_region <- junction_filter(junctions,
    count_thresh = NULL,
    n_samp = NULL,
    types = NULL,
    regions = GRanges(seqnames = "21", ranges = "1-10000000", strand = "*")
)

by_all <- junction_filter(junctions,
    count_thresh = c("raw" = 5),
    n_samp = c("raw" = 1),
    width_range = c(20, 1000),
    types = c("ambig_gene", "unannotated"),
    regions = GRanges(seqnames = "21", ranges = "1-10000000", strand = "*")
)

test_that("junction_filter has correct output", {

    # count
    expect_false(junctions %>%
        assays() %>%
        .[["raw"]] %>%
        apply(
            MARGIN = 1,
            function(x) {
                (sum(x >= 5)) >= 1
            }
        ) %>%
        all())

    expect_true(by_count %>%
        assays() %>%
        .[["raw"]] %>%
        apply(
            MARGIN = 1,
            function(x) {
                (sum(x >= 5)) >= 1
            }
        ) %>%
        all())

    # width
    expect_false(all(width(junctions) <= 1000 & width(junctions) >= 20))
    expect_true(all(width(by_width) <= 1000 & width(by_width) >= 20))

    # type
    expect_true(any(mcols(junctions)[["type"]] %in% c("ambig_gene", "unannotated")))
    expect_false(any(mcols(by_type)[["type"]] %in% c("ambig_gene", "unannotated")))

    # region
    expect_true(length(findOverlaps(junctions, GRanges(seqnames = "21", ranges = "1-10000000", strand = "*"))) > 0)
    expect_true(length(findOverlaps(by_region, GRanges(seqnames = "21", ranges = "1-10000000", strand = "*"))) == 0)

    # all
    expect_true(by_all %>%
        assays() %>%
        .[["raw"]] %>%
        apply(
            MARGIN = 1,
            function(x) {
                (sum(x >= 5)) >= 1
            }
        ) %>%
        all())
    expect_true(all(width(by_all) <= 1000 & width(by_all) >= 20))
    expect_false(any(mcols(by_all)[["type"]] %in% c("ambig_gene", "unannotated")))
    expect_true(length(findOverlaps(by_all, GRanges(seqnames = "21", ranges = "1-10000000", strand = "*"))) == 0)
})
