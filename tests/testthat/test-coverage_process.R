context("Testing coverage processing")

# Load data ---------------------------------------------------------------

# use Genomic state to load txdb (GENCODE v31)
ref <- GenomicState::GenomicStateHub(version = "31", genome = "hg38", filetype = "TxDb")[[1]]

# filter junctions to save time for annotating
# whilst preserving both annotated and unannotated types
# extract random set of 1000 junctions
set.seed(32)

junctions_processed_path <- file.path(tempdir(), "junctions_processed.rda")

if (file.exists(junctions_processed_path)) {
    load(junctions_processed_path)
} else {
    set.seed(32)

    junctions_processed <-
        junction_process(
            junctions_example[sample(seq_len(dim(junctions_example)[1]), 10000), ],
            count_thresh = NULL,
            n_samp = NULL,
            ref,
            sd_const = 0.02
        )

    save(junctions_processed,
        file = junctions_processed_path
    )
}

junctions <- junctions_processed

# Loading and normalising coverage ---------------------------------------------

##### .coverage_exon_intron #####

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

##### .coverage_load_samp #####

# set up a load function avoid having to actually load coverage
# as this is tested itself in test-utils.R for .load_coverage
# uses regions and coverage_path to determine the random seed
# therefore allows testing that structure of output has not been mixed up
# across regions/paths
load_rand <- function(regions, coverage_path, chr_format, sum_fun) {
    seed <- regions %>%
        start() %>%
        mean() %>%
        as.integer()
    seed <- seed + coverage_path

    set.seed(seed)

    cov_rand <- sample(1:100,
        nrow(regions %>% as.data.frame()),
        replace = TRUE
    ) %>%
        as.double()

    return(cov_rand)
}

# set up pseudo paths
# the value of these will be used as random seeds in .load_rand
coverage_paths_case <- c(1, 2)
coverage_paths_control <- c(3, 4, 5)

coverage <- .coverage_load_samp(coverage_regions,
    coverage_paths_case,
    coverage_paths_control,
    coverage_chr_control = "chr",
    load_func = load_rand,
    bp_param = BiocParallel::SerialParam()
)

test_that(".coverage_load_samp output looks correct", {
    expect_true(is(coverage, "list"))
    expect_identical(names(coverage), c("case", "control"))
    expect_identical(
        names(coverage[["case"]]),
        c("exon_coverage_start", "exon_coverage_end", "intron_coverage", "norm_coverage")
    )
    expect_identical(
        names(coverage[["control"]]),
        c("exon_coverage_start", "exon_coverage_end", "intron_coverage", "norm_coverage")
    )
    expect_identical(
        dim(coverage[["case"]][["exon_coverage_start"]]),
        c(length(junctions), length(coverage_paths_case))
    )
    expect_identical(
        dim(coverage[["control"]][["exon_coverage_start"]]),
        c(length(junctions), length(coverage_paths_control))
    )

    expect_true(identical(
        load_rand(
            regions = coverage_regions[["exon_coords_start"]],
            coverage_path = coverage_paths_case[1],
            sum_fun = "mean"
        ),
        coverage[["case"]][["exon_coverage_start"]][, 1]
    ))

    expect_true(identical(
        load_rand(
            regions = coverage_regions[["exon_coords_end"]],
            coverage_path = coverage_paths_case[1],
            sum_fun = "mean"
        ),
        coverage[["case"]][["exon_coverage_end"]][, 1]
    ))

    expect_true(identical(
        load_rand(
            regions = coverage_regions[["intron_coords"]],
            coverage_path = coverage_paths_case[2],
            sum_fun = "mean"
        ),
        coverage[["case"]][["intron_coverage"]][, 2]
    ))

    expect_true(identical(
        load_rand(
            regions = coverage_regions[["norm_coords"]],
            coverage_path = coverage_paths_case[2],
            sum_fun = "sum"
        ),
        coverage[["case"]][["norm_coverage"]][, 2]
    ))
})

##### .coverage_norm #####

coverage_normalised <- .coverage_norm(coverage, norm_const = 2)

test_that(".coverage_norm output looks correct", {
    expect_equivalent(
        coverage_normalised[["case"]][["exon_coverage_start"]][1, ],
        coverage[["case"]][["exon_coverage_start"]][1, ] /
            (coverage[["case"]][["norm_coverage"]][1, ] + 2)
    )

    expect_equivalent(
        coverage_normalised[["case"]][["exon_coverage_start"]][2, ],
        coverage[["case"]][["exon_coverage_start"]][2, ] /
            (coverage[["case"]][["norm_coverage"]][2, ] + 2)
    )
})

##### coverage_norm #####

test_that("coverage_norm catches user-input errors", {
    expect_error(
        coverage_norm(junctions, ref,
            unannot_width,
            coverage_paths_case = c("needs", "to", "be", "length", "2"),
            coverage_paths_control = c("path1", "path2"),
            coverage_chr_control = "chr"
        ),
        "Number of cases must equal the length of coverage_paths_case"
    )

    expect_error(
        coverage_norm(junctions, ref,
            unannot_width,
            coverage_paths_case = c("path1", "path2"),
            coverage_paths_control = c("needs to be at least 2"),
            coverage_chr_control = "chr"
        ),
        "coverage_paths_control must cover at least 2 controls"
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

# Scoring coverage ---------------------------------------------

##### .coverage_score #####

coverage_scores <- .coverage_score(coverage_normalised, .zscore, sd_const = 0.02)

test_that(".coverage_score output generally looks correct", {
    expect_true(is(coverage_scores, "list"))
    expect_false(any(!is.finite(coverage_scores %>% unlist() %>% unlist())))
    expect_identical(dim(coverage_scores[[1]]), dim(coverage_normalised[[1]][[1]]))
})

test_that("zscore has been calculated correctly", {
    control_mean <- coverage_normalised[["control"]][[1]] %>%
        apply(
            MARGIN = 1,
            FUN = mean
        )

    control_sd <- coverage_normalised[["control"]][[1]] %>%
        apply(
            MARGIN = 1,
            FUN = sd
        )

    expect_equivalent(
        coverage_scores[[1]][, 1],
        (coverage_normalised[["case"]][[1]][, 1] - control_mean) / (control_sd + 0.02)
    )

    expect_equivalent(
        coverage_scores[[1]][, 2],
        (coverage_normalised[["case"]][[1]][, 2] - control_mean) / (control_sd + 0.02)
    )
})

##### .coverage_score_max #####

coverage_region_scores_max <- .coverage_score_max(coverage_scores)

coverage_scores_per_samp <- dplyr::tibble(
    exon_coverage_score_start = coverage_scores[["exon_coverage_score_start"]][, 1],
    exon_coverage_score_end = coverage_scores[["exon_coverage_score_end"]][, 1],
    intron_coverage_score = coverage_scores[["intron_coverage_score"]][, 1]
)

test_that(".coverage_score_max output generally looks correct", {
    expect_true(is(coverage_region_scores_max, "list"))
    expect_true(all(unlist(coverage_region_scores_max[["regions"]]) %in% c(1, 2, 3)))
    expect_identical(dim(coverage_region_scores_max[[1]]), dim(coverage_region_scores_max[[2]]))
    expect_identical(dim(coverage_region_scores_max[[1]]), dim(coverage_scores[[1]]))

    expect_identical(
        (abs(coverage_scores_per_samp) == abs(coverage_region_scores_max[["scores"]][, 1])) %>%
            apply(
                MARGIN = 1,
                FUN = which
            ) %>%
            lapply(FUN = min) %>% # takes the min index if multiple regions have same score
            unlist() %>%
            unname(),
        coverage_region_scores_max[["regions"]][, 1]
    )
})

##### coverage_score #####

names(coverage[["case"]]) <- c("not", "correct", "names")

test_that("coverage_score catches user-input errors", {
    expect_error(
        coverage_score(coverage = coverage[[1]]),
        "coverage should have the names 'case' and 'control'"
    )

    expect_error(
        coverage_score(coverage = coverage),
        "coverage matrices should be named 'exon_coverage_start', 'exon_coverage_end', 'intron_coverage'"
    )
})
