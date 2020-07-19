context("Testing normalistion of junction counts")

##### junction_norm #####

junctions <- junction_norm(junctions_example)

test_that("junction_norm general output looks correct", {
    expect_true(methods::isClass(junctions, "RangedSummarizedExperiment"))
    expect_equal(
        colnames(SummarizedExperiment::assays(junctions)[["raw"]]),
        colnames(SummarizedExperiment::assays(junctions)[["norm"]])
    )
    expect_equal(
        nrow(SummarizedExperiment::assays(junctions)[["raw"]]),
        nrow(SummarizedExperiment::assays(junctions)[["norm"]])
    )
    expect_false(any(lengths(mcols(junctions)[["clusters"]]) == 0))
    expect_identical(
        names(mcols(junctions)[["clusters"]]),
        mcols(junctions)[["index"]] %>% as.character()
    )
})

# all clusters with only a single junction should have a value of 1 or 0
mono_cluster <- which(lengths(mcols(junctions)[["clusters"]]) == 1)
mono_cluster_counts <- SummarizedExperiment::assays(junctions)[["norm"]][mono_cluster, ]

test_that("normalised counts have expected values", {
    expect_false(any(is.na(SummarizedExperiment::assays(junctions)[["norm"]])))
    expect_true(all(SummarizedExperiment::assays(junctions)[["norm"]] <= 1))
    expect_true(all(SummarizedExperiment::assays(junctions)[["norm"]] >= 0))
    expect_identical(
        which(SummarizedExperiment::assays(junctions)[["raw"]] == 0),
        which(SummarizedExperiment::assays(junctions)[["norm"]] == 0)
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
            expect_cluster_sum_counts <- SummarizedExperiment::assays(junctions)[["raw"]][expect_cluster, ] %>%
                apply(MARGIN = 2, FUN = sum)
        } else {
            expect_cluster_sum_counts <- SummarizedExperiment::assays(junctions)[["raw"]][expect_cluster, ]
        }

        expect_norm_counts <- SummarizedExperiment::assays(junction_start_end_to_test[["start"]])[["raw"]] / expect_cluster_sum_counts
        expect_norm_counts[is.na(expect_norm_counts)] <- 0

        check <- all(check, identical(
            expect_cluster,
            mcols(junctions)[["clusters"]][[i]]
        ))

        check <- all(check, identical(
            expect_norm_counts[1, ],
            SummarizedExperiment::assays(junctions)[["norm"]][i, ]
        ))
    }

    return(check)
}

test_that("raw counts have been correctly normalised", {
    expect_true(norm_check(junctions, 50))
})
