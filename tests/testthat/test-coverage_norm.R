context("Testing the loading and normalisation of coverage")

github <- TRUE

##### .coverage_exon_intron #####

junctions <- junctions_annot_example[1:1000]

ref <- "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"

suppressWarnings(expr = {
    ref <- GenomicFeatures::makeTxDbFromGFF(ref)
})

coverage_regions <- .coverage_exon_intron(junctions, unannot_width = 10)

# only testing annotated/unannotated
# as .coverage_region operates equally across all types
# and these two are easiest to test
annot_indexes <- which(mcols(junctions)[["type"]] == "annotated")
unannot_indexes <- which(mcols(junctions)[["type"]] == "unannotated")

exon_widths_start <- mcols(junctions)[["exon_width_start"]][annot_indexes] %>%
    min() %>%
    unname()

exon_widths_end <- mcols(junctions)[["exon_width_end"]][annot_indexes] %>%
    min() %>%
    unname()

test_that(".coverage_exon_intron output looks correct", {
    expect_equivalent(
        end(coverage_regions[["exon_coords_start"]]),
        start(coverage_regions[["intron_coords"]]) - 1
    )
    expect_equivalent(
        start(coverage_regions[["exon_coords_end"]]),
        end(coverage_regions[["intron_coords"]]) + 1
    )

    expect_equivalent(
        (end(coverage_regions[["exon_coords_start"]][annot_indexes]) - (exon_widths_start - 1)),
        start(coverage_regions[["exon_coords_start"]][annot_indexes])
    )
    expect_equivalent(
        (start(coverage_regions[["exon_coords_end"]][annot_indexes]) + (exon_widths_end - 1)),
        end(coverage_regions[["exon_coords_end"]][annot_indexes])
    )

    expect_true(all(width(coverage_regions[["exon_coords_start"]][unannot_indexes]) == 10))
    expect_true(all(width(coverage_regions[["exon_coords_end"]][unannot_indexes]) == 10))
})

##### .coverage_norm_region #####

coverage_regions <- .coverage_norm_region(junctions, ref, coverage_regions)

ref_genes <- ref %>%
    GenomicFeatures::genes(columns = c("gene_id", "exon_name"))

coverage_norm_region_check <- function(junctions, coverage_regions, ref_genes, n) {
    check <- TRUE

    for (i in sample(seq_along(junctions), n)) {
        junction_to_test <- junctions[i]

        gene_ids <- mcols(junction_to_test)[["gene_id_junction"]] %>%
            unlist()

        if (length(gene_ids) == 0) {
            check <- all(check, identical(start(coverage_regions[["norm_coords"]][i]), start(coverage_regions[["exon_coords_start"]][i])))
            check <- all(check, identical(end(coverage_regions[["norm_coords"]][i]), end(coverage_regions[["exon_coords_end"]][i])))
            check <- all(check, identical(
                seqnames(coverage_regions[["norm_coords"]][i]) %>% as.character(),
                seqnames(junction_to_test) %>% as.character()
            ))
            check <- all(check, identical(strand(coverage_regions[["norm_coords"]][i]), strand(junction_to_test)))
        } else {
            gene_to_test <- ref_genes[mcols(ref_genes)[["gene_id"]] %in% gene_ids]
            gene_to_test <- gene_to_test[which.min(width(gene_to_test))]

            check <- all(check, length(findOverlaps(gene_to_test, coverage_regions[["norm_coords"]][i], type = "equal")) == 1)
        }

        if (check == FALSE) {
            stop(i)
        }
    }

    return(check)
}

test_that(".coverage_norm_region output looks correct", {
    expect_true(coverage_norm_region_check(junctions, coverage_regions, ref_genes, n = 50))
})

##### .coverage_case_control_load #####

coverage_paths_case <- list.files("/data/RNA_seq_diag/mito/bw/", full.names = T)[1:2]
coverage_paths_control <- list.files("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/", full.names = T)[1]

test_that(".coverage_case_control_load output looks correct", {
    if (github) {
        skip("skipping as not testing loading coverage from remote files yet")
    }

    case_control_coverage <- .coverage_case_control_load(coverage_regions,
        coverage_paths_case,
        coverage_paths_control,
        coverage_chr_control = "chr"
    )

    expect_true(is(case_control_coverage, "list"))
    expect_identical(names(case_control_coverage), c("case", "control"))
    expect_identical(
        names(case_control_coverage[["case"]]),
        c("exon_coverage_start", "exon_coverage_end", "intron_coverage", "norm_coverage")
    )
    expect_identical(
        names(case_control_coverage[["control"]]),
        c("exon_coverage_start", "exon_coverage_end", "intron_coverage", "norm_coverage")
    )
    expect_identical(
        dim(case_control_coverage[["case"]][["exon_coverage_start"]]),
        c(length(junctions), length(coverage_paths_case))
    )
    expect_identical(
        dim(case_control_coverage[["control"]][["exon_coverage_start"]]),
        c(length(junctions), length(coverage_paths_control))
    )

    expect_true(identical(
        .coverage_load(coverage_paths_case[1],
            regions = coverage_regions[["exon_coords_start"]],
            sum_fun = "mean"
        ),
        case_control_coverage[["case"]][["exon_coverage_start"]][, 1]
    ))

    expect_true(identical(
        .coverage_load(coverage_paths_case[1],
            regions = coverage_regions[["exon_coords_end"]],
            sum_fun = "mean"
        ),
        case_control_coverage[["case"]][["exon_coverage_end"]][, 1]
    ))

    expect_true(identical(
        .coverage_load(coverage_paths_case[2],
            regions = coverage_regions[["intron_coords"]],
            sum_fun = "mean"
        ),
        case_control_coverage[["case"]][["intron_coverage"]][, 2]
    ))

    expect_true(identical(
        .coverage_load(coverage_paths_case[2],
            regions = coverage_regions[["norm_coords"]],
            sum_fun = "sum"
        ),
        case_control_coverage[["case"]][["norm_coverage"]][, 2]
    ))
})

##### .coverage_norm #####

test_that(".coverage_norm output looks correct", {
    if (github) {
        skip("skipping as not testing loading coverage from remote files yet")
    }

    # to be removed (as repeat of above)
    # when can be taken out into global env
    # after remote coverage loading has been tested
    case_control_coverage <- .coverage_case_control_load(coverage_regions,
        coverage_paths_case,
        coverage_paths_control,
        coverage_chr_control = "chr"
    )

    case_control_coverage_norm <- .coverage_norm(case_control_coverage)

    # no non-finite such as Inf caused by 0 denominator
    expect_false(any(!is.finite(case_control_coverage_norm %>% unlist() %>% unlist())))

    expect_identical(
        case_control_coverage_norm[["case"]][["exon_coverage_start"]][1, ],
        case_control_coverage[["case"]][["exon_coverage_start"]][1, ] / case_control_coverage[["case"]][["norm_coverage"]][1, ]
    )

    expect_identical(
        case_control_coverage_norm[["case"]][["exon_coverage_start"]][2, ],
        case_control_coverage[["case"]][["exon_coverage_start"]][2, ] / case_control_coverage[["case"]][["norm_coverage"]][2, ]
    )
})

##### coverage_norm #####

test_that("coverage_norm catches user-input errors", {
    expect_error(
        coverage_norm(junctions, ref,
            unannot_width,
            coverage_paths_case = c("needs", "to", "be", "2"),
            coverage_paths_control = c("path1", "path2"),
            coverage_chr_control = "chr"
        ),
        "Number of cases must equal the length of coverage_paths_case"
    )

    expect_error(
        coverage_norm(junctions, ref,
            unannot_width,
            coverage_paths_case = c("path1", "path2"),
            coverage_paths_control = c("path"),
            coverage_chr_control = "chr"
        ),
        "coverage_paths_control must coverageer at least 2 controls"
    )

    expect_error(
        coverage_norm(junctions, ref,
            unannot_width,
            coverage_paths_case = c("path1", "path2"),
            coverage_paths_control = c("path1", "path2"),
            coverage_chr_control = "not_a_chr"
        ),
        "coverage_chr_control must be one of 'chr' or 'no_chr'"
    )
})
