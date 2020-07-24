context("Testing the loading and normalisation of coverage")

##### .cov_exon_intron #####

junctions <- junctions_annot_example[1:1000]

ref <- "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"

suppressWarnings(expr = {
    ref <- GenomicFeatures::makeTxDbFromGFF(ref)
})

cov_regions <- .cov_exon_intron(junctions, unannot_width = 10)

# only testing annotated/unannotated
# as .cov_region operates equally across all types
# and these two are easiest to test
annot_indexes <- which(mcols(junctions)[["type"]] == "annotated")
unannot_indexes <- which(mcols(junctions)[["type"]] == "unannotated")

exon_widths_start <- mcols(junctions)[["exon_width_start"]][annot_indexes] %>%
    min() %>%
    unname()

exon_widths_end <- mcols(junctions)[["exon_width_end"]][annot_indexes] %>%
    min() %>%
    unname()

test_that(".cov_exon_intron output looks correct", {
    expect_equivalent(
        end(cov_regions[["exon_coords_start"]]),
        start(cov_regions[["intron_coords"]]) - 1
    )
    expect_equivalent(
        start(cov_regions[["exon_coords_end"]]),
        end(cov_regions[["intron_coords"]]) + 1
    )

    expect_equivalent(
        (end(cov_regions[["exon_coords_start"]][annot_indexes]) - (exon_widths_start - 1)),
        start(cov_regions[["exon_coords_start"]][annot_indexes])
    )
    expect_equivalent(
        (start(cov_regions[["exon_coords_end"]][annot_indexes]) + (exon_widths_end - 1)),
        end(cov_regions[["exon_coords_end"]][annot_indexes])
    )

    expect_true(all(width(cov_regions[["exon_coords_start"]][unannot_indexes]) == 10))
    expect_true(all(width(cov_regions[["exon_coords_end"]][unannot_indexes]) == 10))
})

##### .cov_norm_region #####

cov_regions <- .cov_norm_region(junctions, ref, cov_regions)

ref_genes <- ref %>%
    GenomicFeatures::genes(columns = c("gene_id", "exon_name"))

cov_norm_region_check <- function(junctions, cov_regions, ref_genes, n) {
    check <- TRUE

    for (i in sample(seq_along(junctions), n)) {
        junction_to_test <- junctions[i]

        gene_ids <- mcols(junction_to_test)[["gene_id_junction"]] %>%
            unlist()

        if (length(gene_ids) == 0) {
            check <- all(check, identical(start(cov_regions[["norm_coords"]][i]), start(cov_regions[["exon_coords_start"]][i])))
            check <- all(check, identical(end(cov_regions[["norm_coords"]][i]), end(cov_regions[["exon_coords_end"]][i])))
            check <- all(check, identical(
                seqnames(cov_regions[["norm_coords"]][i]) %>% as.character(),
                seqnames(junction_to_test) %>% as.character()
            ))
            check <- all(check, identical(strand(cov_regions[["norm_coords"]][i]), strand(junction_to_test)))
        } else {
            gene_to_test <- ref_genes[mcols(ref_genes)[["gene_id"]] %in% gene_ids]
            gene_to_test <- gene_to_test[which.min(width(gene_to_test))]

            check <- all(check, length(findOverlaps(gene_to_test, cov_regions[["norm_coords"]][i], type = "equal")) == 1)
        }

        if (check == FALSE) {
            stop(i)
        }
    }

    return(check)
}

test_that(".cov_norm_region output looks correct", {
    expect_true(cov_norm_region_check(junctions, cov_regions, ref_genes, n = 50))
})

##### .cov_case_control_load #####

cov_paths_case <- list.files("/data/RNA_seq_diag/mito/bw/", full.names = T)[1:2]
cov_paths_control <- list.files("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/", full.names = T)[1]

case_control_cov <- .cov_case_control_load(cov_regions,
    cov_paths_case,
    cov_paths_control,
    cov_chr_control = "chr"
)

github <- TRUE

test_that(".cov_case_control_load output looks correct", {
    expect_true(is(case_control_cov, "list"))
    expect_identical(names(case_control_cov), c("case", "control"))
    expect_identical(
        names(case_control_cov[["case"]]),
        c("exon_cov_start", "exon_cov_end", "intron_cov", "norm_cov")
    )
    expect_identical(
        names(case_control_cov[["control"]]),
        c("exon_cov_start", "exon_cov_end", "intron_cov", "norm_cov")
    )
    expect_identical(
        dim(case_control_cov[["case"]][["exon_cov_start"]]),
        c(length(junctions), length(cov_paths_case))
    )
    expect_identical(
        dim(case_control_cov[["control"]][["exon_cov_start"]]),
        c(length(junctions), length(cov_paths_control))
    )

    if (github) {
        skip()
    }

    expect_true(identical(
        .cov_load(cov_paths_case[1],
            regions = cov_regions[["exon_coords_start"]],
            sum_fun = "mean"
        ),
        case_control_cov[["case"]][["exon_cov_start"]][, 1]
    ))

    expect_true(identical(
        .cov_load(cov_paths_case[1],
            regions = cov_regions[["exon_coords_end"]],
            sum_fun = "mean"
        ),
        case_control_cov[["case"]][["exon_cov_end"]][, 1]
    ))

    expect_true(identical(
        .cov_load(cov_paths_case[2],
            regions = cov_regions[["intron_coords"]],
            sum_fun = "mean"
        ),
        case_control_cov[["case"]][["intron_cov"]][, 2]
    ))

    expect_true(identical(
        .cov_load(cov_paths_case[2],
            regions = cov_regions[["norm_coords"]],
            sum_fun = "sum"
        ),
        case_control_cov[["case"]][["norm_cov"]][, 2]
    ))
})

##### .cov_norm #####

case_control_cov_norm <- .cov_norm(case_control_cov)

test_that(".cov_norm output looks correct", {
    expect_identical(
        case_control_cov_norm[["case"]][["exon_cov_start"]][1, ],
        case_control_cov[["case"]][["exon_cov_start"]][1, ] / case_control_cov[["case"]][["norm_cov"]][1, ]
    )

    expect_identical(
        case_control_cov_norm[["case"]][["exon_cov_start"]][2, ],
        case_control_cov[["case"]][["exon_cov_start"]][2, ] / case_control_cov[["case"]][["norm_cov"]][2, ]
    )
})

##### junction_cov_norm #####

test_that("junction_cov_norm catches user-input errors", {
    expect_error(
        junction_cov_norm(junctions, ref,
            unannot_width,
            cov_paths_case = c("needs", "to", "be", "2"),
            cov_paths_control = c("path"),
            cov_chr_control = "chr"
        ),
        "Number of cases must equal the length of cov_paths_case"
    )

    expect_error(
        junction_cov_norm(junctions, ref,
            unannot_width,
            cov_paths_case,
            cov_paths_control,
            cov_chr_control = "not_a_chr"
        ),
        "cov_chr_control must be one of 'chr' or 'no_chr'"
    )
})
