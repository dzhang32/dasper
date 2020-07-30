context("Test functions that score coverage")

github <- TRUE

##### .coverage_score #####

if (!github) {
    ref <- "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
    suppressWarnings(expr = {
        ref <- GenomicFeatures::makeTxDbFromGFF(ref)
    })
    junctions <- junction_norm(junctions_annot_example)
    junctions <- junction_score(junctions)
    coverage_paths_case <- list.files("/data/RNA_seq_diag/mito/bw/", full.names = T)[1:2]
    coverage_paths_control <- list.files("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/", full.names = T)[1:2]
    coverage <- coverage_norm(junctions, ref, unannot_width = 20, coverage_paths_case, coverage_paths_control, coverage_chr_control = "chr")

    coverage_scores <- .coverage_score(coverage, .zscore, sd_const = 0.02)
}

test_that(".coverage_score output generally looks correct", {
    if (github) {
        skip("skipping as not testing loading coverage from remote files yet")
    }

    expect_true(is(coverage_scores, "list"))
    expect_false(any(!is.finite(coverage_scores %>% unlist() %>% unlist())))
    expect_identical(dim(coverage_scores[[1]]), dim(coverage[[1]][[1]]))
})

test_that("zscore has been calculated correctly", {
    if (github) {
        skip("skipping as not testing loading coverage from remote files yet")
    }

    control_mean <- coverage[["control"]][[1]] %>%
        apply(
            MARGIN = 1,
            FUN = mean
        )

    control_sd <- coverage[["control"]][[1]] %>%
        apply(
            MARGIN = 1,
            FUN = sd
        )

    expect_equivalent(
        coverage_scores[[1]][, 1],
        (coverage[["case"]][[1]][, 1] - control_mean) / (control_sd + 0.02)
    )

    expect_equivalent(
        coverage_scores[[1]][, 2],
        (coverage[["case"]][[1]][, 2] - control_mean) / (control_sd + 0.02)
    )
})

##### .coverage_score_max #####

test_that(".coverage_score_max output generally looks correct", {
    if (github) {
        skip("skipping as not testing loading coverage from remote files yet")
    }

    coverage_region_scores_max <- .coverage_score_max(coverage_scores)

    coverage_scores_per_samp <- dplyr::tibble(
        exon_coverage_score_start = coverage_scores[["exon_coverage_score_start"]][, 1],
        exon_coverage_score_end = coverage_scores[["exon_coverage_score_end"]][, 1],
        intron_coverage_score = coverage_scores[["intron_coverage_score"]][, 1]
    )

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
