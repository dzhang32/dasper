context("Test functions that score coverage")

github <- TRUE

##### .cov_score #####

if (!github) {
    ref <- "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
    suppressWarnings(expr = {
        ref <- GenomicFeatures::makeTxDbFromGFF(ref)
    })
    junctions <- junction_norm(junctions_annot_example)
    junctions <- junction_score(junctions)
    cov_paths_case <- list.files("/data/RNA_seq_diag/mito/bw/", full.names = T)[1:2]
    cov_paths_control <- list.files("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/", full.names = T)[1:2]
    cov <- junction_cov_norm(junctions, ref, unannot_width = 20, cov_paths_case, cov_paths_control, cov_chr_control = "chr")

    cov_scores <- .cov_score(cov, .zscore, sd_const = 0.02)
}

test_that(".cov_score output generally looks correct", {
    if (github) {
        skip("skipping as not testing loading coverage from remote files yet")
    }

    expect_true(is(cov_scores, "list"))
    expect_false(any(!is.finite(cov_scores %>% unlist() %>% unlist())))
    expect_identical(dim(cov_scores[[1]]), dim(cov[[1]][[1]]))
})

test_that("zscore has been calculated correctly", {
    if (github) {
        skip("skipping as not testing loading coverage from remote files yet")
    }

    control_mean <- cov[["control"]][[1]] %>%
        apply(
            MARGIN = 1,
            FUN = mean
        )

    control_sd <- cov[["control"]][[1]] %>%
        apply(
            MARGIN = 1,
            FUN = sd
        )

    expect_equivalent(
        cov_scores[[1]][, 1],
        (cov[["case"]][[1]][, 1] - control_mean) / (control_sd + 0.02)
    )

    expect_equivalent(
        cov_scores[[1]][, 2],
        (cov[["case"]][[1]][, 2] - control_mean) / (control_sd + 0.02)
    )
})

##### .cov_score_max #####

test_that(".cov_score_max output generally looks correct", {
    if (github) {
        skip("skipping as not testing loading coverage from remote files yet")
    }

    cov_region_scores_max <- .cov_score_max(cov_scores)

    cov_scores_per_samp <- dplyr::tibble(
        exon_cov_score_start = cov_scores[["exon_cov_score_start"]][, 1],
        exon_cov_score_end = cov_scores[["exon_cov_score_end"]][, 1],
        intron_cov_score = cov_scores[["intron_cov_score"]][, 1]
    )

    expect_true(is(cov_region_scores_max, "list"))
    expect_true(all(unlist(cov_region_scores_max[["regions"]]) %in% c(1, 2, 3)))
    expect_identical(dim(cov_region_scores_max[[1]]), dim(cov_region_scores_max[[2]]))
    expect_identical(dim(cov_region_scores_max[[1]]), dim(cov_scores[[1]]))

    expect_identical(
        (abs(cov_scores_per_samp) == abs(cov_region_scores_max[["scores"]][, 1])) %>%
            apply(
                MARGIN = 1,
                FUN = which
            ) %>%
            lapply(FUN = min) %>% # takes the min index if multiple regions have same score
            unlist() %>%
            unname(),
        cov_region_scores_max[["regions"]][, 1]
    )
})

##### junction_cov_score #####

names(cov[["case"]]) <- c("not", "correct", "names")

test_that("junction_cov_score catches user-input errors", {
    expect_error(
        junction_cov_score(cov = cov[[1]]),
        "cov should have the names 'case' and 'control'"
    )

    expect_error(
        junction_cov_score(cov = cov),
        "coverage matrices should be named 'exon_cov_start', 'exon_cov_end', 'intron_cov'"
    )
})
