context("Testing junction processing")

# Annotating junctions ----------------------------------------------------

##### .junction_annot_tidy #####

x <- GRanges(c("1:1-100", "2:2-200", "3:3-300", "4:4-400"))
strand(x) <- c("+", "-", "*", "*")
mcols(x)["strand_start"] <- CharacterList("1" = "+", "2" = "*", "3" = "+", "4" = "*")
mcols(x)["strand_end"] <- mcols(x)["strand_start"]
strand_tidy <- as.character(strand(.junction_annot_tidy(x, cols_to_merge = "strand")))

test_that(".junction_annot_tidy correctly infers strand", {
    expect_identical(strand_tidy, c("+", "-", "+", "*"))
})

##### junction_annot #####

ref <- "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
suppressWarnings(expr = {
    ref <- .ref_load(ref)
})

junctions <- junction_annot(junctions_example, ref)

test_that("junction_annot output generally looks correct", {
    expect_true(methods::isClass(junctions, "RangedSummarisedExperiment"))
    expect_false(any(is.na(mcols(junctions)[["type"]])))
    expect_identical(
        lengths(mcols(junctions)[["exon_name_start"]]),
        lengths(mcols(junctions)[["exon_width_start"]])
    )
    expect_identical(
        lengths(mcols(junctions)[["exon_name_end"]]),
        lengths(mcols(junctions)[["exon_width_end"]])
    )
})

test_that("junction_annot catches user-input errors", {
    expect_error(
        junction_annot(junctions = "not_a_RSE", ref = ref),
        "junctions must be in a RangedSummarisedExperiment format"
    )
    expect_error(
        junction_annot(junctions = junctions, ref = 50),
        "ref must either be a path to the .gtf/gff3 file or a pre-loaded TxDb object"
    )
})

# split junctions by categories
annotated <- junctions[mcols(junctions)[["type"]] == "annotated"]
novel_donor <- junctions[mcols(junctions)[["type"]] == "novel_donor"]
novel_acceptor <- junctions[mcols(junctions)[["type"]] == "novel_acceptor"]
novel_exon_skip <- junctions[mcols(junctions)[["type"]] == "novel_exon_skip"]
novel_combo <- junctions[mcols(junctions)[["type"]] == "novel_combo"]
ambig_gene <- junctions[mcols(junctions)[["type"]] == "ambig_gene"]
unannotated <- junctions[mcols(junctions)[["type"]] == "unannotated"]

test_that("junction categories meet expectations", {
    expect_identical(
        levels(mcols(junctions)[["type"]]),
        c(
            "annotated",
            "novel_acceptor",
            "novel_donor",
            "novel_exon_skip",
            "novel_combo",
            "ambig_gene",
            "unannotated"
        )
    )

    # annotated
    expect_true(all(lengths(mcols(annotated)[["gene_id_junction"]]) > 0))
    expect_true(all(lengths(mcols(annotated)[["gene_id_start"]]) > 0 &
        lengths(mcols(annotated)[["gene_id_end"]]) > 0))

    # novel_donor
    expect_true(all(lengths(mcols(novel_donor[strand(novel_donor) == "+"])[["gene_id_start"]]) == 0))
    expect_true(all(lengths(mcols(novel_donor[strand(novel_donor) == "+"])[["gene_id_end"]]) > 0))
    expect_true(all(lengths(mcols(novel_donor[strand(novel_donor) == "-"])[["gene_id_start"]]) > 0))
    expect_true(all(lengths(mcols(novel_donor[strand(novel_donor) == "-"])[["gene_id_end"]]) == 0))

    # novel_acceptor
    expect_true(all(lengths(mcols(novel_acceptor[strand(novel_acceptor) == "+"])[["gene_id_start"]]) > 0))
    expect_true(all(lengths(mcols(novel_acceptor[strand(novel_acceptor) == "+"])[["gene_id_end"]]) == 0))
    expect_true(all(lengths(mcols(novel_acceptor[strand(novel_acceptor) == "-"])[["gene_id_start"]]) == 0))
    expect_true(all(lengths(mcols(novel_acceptor[strand(novel_acceptor) == "-"])[["gene_id_end"]]) > 0))

    # novel_combo
    expect_false(any(any(mcols(novel_combo)[["tx_name_start"]] %in% mcols(novel_combo)[["tx_name_end"]])))

    # novel_exon_skip
    expect_true(all(any(mcols(novel_exon_skip)[["tx_name_start"]] %in% mcols(novel_exon_skip)[["tx_name_end"]])))

    # ambig_gene
    expect_true(all(lengths(ambig_gene$gene_id_junc) > 1))

    # unannotated
    expect_false(any(mcols(unannotated)[["in_ref"]]))
    expect_true(length(c(
        unlist(mcols(unannotated)[["exon_name_start"]]),
        unlist(mcols(unannotated)[["exon_name_end"]])
    )) == 0)
})

ref_exons <- ref %>% GenomicFeatures::exons(columns = c("gene_id", "tx_name", "exon_name"))

# define function to manually check n randomly sampled junctions
# to make exon annotation from junction_annot matches expections
annot_check <- function(junctions, ref_exons, n) {
    check <- TRUE

    start(junctions) <- start(junctions) - 1
    end(junctions) <- end(junctions) + 1

    junctions_start_end <- .get_start_end(junctions)
    ref_exons_start_end <- .get_start_end(ref_exons)

    for (i in sample(1:length(junctions), n, replace = TRUE)) {
        junction_to_test <- junctions_start_end %>%
            lapply(FUN = function(x) {
                x[i]
            })

        expect_exons_start <-
            ref_exons[findOverlaps(ref_exons_start_end[["end"]],
                junction_to_test[["start"]],
                ignore.strand = FALSE
            ) %>%
                queryHits()]

        expect_exons_end <-
            ref_exons[findOverlaps(ref_exons_start_end[["start"]],
                junction_to_test[["end"]],
                ignore.strand = FALSE
            ) %>%
                queryHits()]

        # check exons
        check <- all(check, identical(
            expect_exons_start$exon_name %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["exon_name_start"]] %>% unlist() %>% unname() %>% sort()
        ))

        check <- all(check, identical(
            expect_exons_end$exon_name %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["exon_name_end"]] %>% unlist() %>% unname() %>% sort()
        ))

        # check exon widths
        check <- all(check, identical(
            expect_exons_start %>% width() %>% sort(),
            mcols(junction_to_test[["start"]])[["exon_width_start"]] %>% unlist() %>% unname() %>% sort()
        ))

        check <- all(check, identical(
            expect_exons_end %>% width() %>% sort(),
            mcols(junction_to_test[["start"]])[["exon_width_end"]] %>% unlist() %>% unname() %>% sort()
        ))

        # check transcripts
        check <- all(check, identical(
            expect_exons_start$tx_name %>% unlist() %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["tx_name_start"]] %>% unlist() %>% unname() %>% sort()
        ))

        check <- all(check, identical(
            expect_exons_end$tx_name %>% unlist() %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["tx_name_end"]] %>% unlist() %>% unname() %>% sort()
        ))

        # check genes
        check <- all(check, identical(
            expect_exons_start$gene_id %>% unlist() %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["gene_id_start"]] %>% unlist() %>% unname() %>% sort()
        ))

        check <- all(check, identical(
            expect_exons_end$gene_id %>% unlist() %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["gene_id_end"]] %>% unlist() %>% unname() %>% sort()
        ))

        if (check == FALSE) {
            stop(print(i))
        }
    }

    return(check)
}

test_that("exon annotation has been correctly retreived", {
    expect_true(annot_check(annotated, ref_exons, 10))
    expect_true(annot_check(novel_donor, ref_exons, 10))
    expect_true(annot_check(novel_acceptor, ref_exons, 10))
    expect_true(annot_check(novel_exon_skip, ref_exons, 10))
    expect_true(annot_check(novel_combo, ref_exons, 5))
    expect_true(annot_check(ambig_gene, ref_exons, 5))
    expect_true(annot_check(unannotated, ref_exons, 10))
})


# Filter junctions --------------------------------------------------------

##### junction_filter #####

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

# Normalise junctions -----------------------------------------------------

##### junction_norm #####

junctions <- junction_norm(junctions)

test_that("junction_norm general output looks correct", {
    expect_true(methods::isClass(junctions, "RangedSummarizedExperiment"))
    expect_equal(
        colnames(assays(junctions)[["raw"]]),
        colnames(assays(junctions)[["norm"]])
    )
    expect_equal(
        nrow(assays(junctions)[["raw"]]),
        nrow(assays(junctions)[["norm"]])
    )
    expect_false(any(lengths(mcols(junctions)[["clusters"]]) == 0))
    expect_identical(
        names(mcols(junctions)[["clusters"]]),
        mcols(junctions)[["index"]] %>% as.character()
    )
})

# all clusters with only a single junction should have a value of 1 or 0
mono_cluster <- which(lengths(mcols(junctions)[["clusters"]]) == 1)
mono_cluster_counts <- assays(junctions)[["norm"]][mono_cluster, ]

test_that("normalised counts have expected values", {
    expect_false(any(is.na(assays(junctions)[["norm"]])))
    expect_true(all(assays(junctions)[["norm"]] <= 1))
    expect_true(all(assays(junctions)[["norm"]] >= 0))
    expect_identical(
        which(assays(junctions)[["raw"]] == 0),
        which(assays(junctions)[["norm"]] == 0)
    )
    expect_true(all(unlist(mono_cluster_counts) %in% c(0, 1)))
})

norm_check <- function(junctions, n) {
    check <- TRUE

    junctions_start_end <- .get_start_end(junctions)

    # only check clusters that contain at least 2 junctions
    polt_cluster <- which(lengths(mcols(junctions)[["clusters"]]) > 1)

    for (i in sample(polt_cluster, n)) {
        junction_start_end_to_test <- junctions_start_end %>%
            lapply(FUN = function(x) {
                x[i]
            })

        start_hits <- findOverlaps(
            junctions_start_end[["start"]],
            junction_start_end_to_test[["start"]]
        )

        end_hits <- findOverlaps(
            junctions_start_end[["end"]],
            junction_start_end_to_test[["end"]]
        )

        expect_cluster <- c(queryHits(start_hits), queryHits(end_hits)) %>%
            unique() %>%
            sort()

        # if the cluster only has 1 junction, we don't need to sum
        if (length(expect_cluster) >= 2) {
            expect_cluster_sum_counts <- assays(junctions)[["raw"]][expect_cluster, ] %>%
                apply(MARGIN = 2, FUN = sum)
        } else {
            expect_cluster_sum_counts <- assays(junctions)[["raw"]][expect_cluster, ]
        }

        expect_norm_counts <- assays(junction_start_end_to_test[["start"]])[["raw"]] / expect_cluster_sum_counts
        expect_norm_counts[is.na(expect_norm_counts)] <- 0

        check <- all(check, identical(
            expect_cluster,
            mcols(junctions)[["clusters"]][[i]]
        ))

        check <- all(check, identical(
            expect_norm_counts[1, ],
            assays(junctions)[["norm"]][i, ]
        ))
    }

    return(check)
}

test_that("raw counts have been correctly normalised", {
    expect_true(norm_check(junctions, 50))
})

# Score junctions ---------------------------------------------------------

##### junction_score #####

junctions_w_score <- junction_score(junctions,
    score_func = .zscore,
    sd_const = 0.02
) # try an sd_const other than default

test_that("junction_score has correct output", {
    expect_true(is(junctions_w_score, "RangedSummarizedExperiment"))
    expect_identical(dim(junctions_w_score)[1], dim(junctions)[1])
    expect_identical(dim(junctions_w_score)[2], sum(colData(junctions)[["case_control"]] == "case"))
    expect_identical(names(assays(junctions_w_score)), c("raw", "norm", "direction", "score"))
    expect_false(any(is.na(assays(junctions_w_score)[["score"]])))
    expect_false(any(!is.finite(assays(junctions_w_score)[["score"]])))
    expect_true(all(unlist(assays(junctions_w_score)[["direction"]]) %in% c(-1, 1)))
})

junctions_no_colData <- junctions_example
SummarizedExperiment::colData(junctions_no_colData)[["case_control"]] <- NULL

test_that("junction_score catches user-input errors", {
    expect_error(
        junction_score(junctions_no_colData),
        "must include the column 'case_control'"
    )

    expect_error(
        junction_score(junctions_example),
        "Junctions must include the 'norm' assay"
    )
})

# vectorised method of performing z-score
# test whether this produces equivalent output to junction_score
case_count <- assays(junctions)[["norm"]][, colData(junctions)[["case_control"]] == "case"]
control_count <- assays(junctions)[["norm"]][, colData(junctions)[["case_control"]] == "control"]

control_mean <-
    control_count %>%
    apply(MARGIN = 1, FUN = mean)

control_sd <-
    control_count %>%
    apply(MARGIN = 1, FUN = sd)

case_score <- (case_count - control_mean) / (control_sd + 0.02)

test_that("direction and score have been calculated correctly", {
    expect_identical(
        assays(junctions_w_score)[["direction"]],
        ifelse(case_score > 0, 1, -1)
    )

    expect_equivalent(
        assays(junctions_w_score)[["score"]],
        case_score %>% unlist()
    )
})

# Process junctions -------------------------------------------------------

##### junction_process #####

# double check all functions fit easily with the pipe
# also that the order of junction_filter/junction_annot
# does not impact result (diff order to junction_process)
junctions_processed <- junctions_example %>%
    junction_annot(ref) %>%
    junction_filter(
        count_thresh = c("raw" = 5),
        n_samp = c("raw" = 1),
        width_range = c(25, 1000000),
        types = c("ambig_gene", "unannotated"),
        regions = GRanges(seqnames = "21", ranges = "1-10000000", strand = "*")
    ) %>%
    junction_norm() %>%
    junction_score(sd_const = 0.02)

junctions_processed_2 <-
    junction_process(
        junctions_example,
        ref,
        count_thresh = c("raw" = 5),
        n_samp = c("raw" = 1),
        width_range = c(25, 1000000),
        types = c("ambig_gene", "unannotated"),
        regions = GRanges(seqnames = "21", ranges = "1-10000000", strand = "*"),
        score_func = .zscore,
        sd_const = 0.02
    )

test_that("junction_process has the correct output", {
    expect_equivalent(junctions_processed, junctions_processed_2)
})
