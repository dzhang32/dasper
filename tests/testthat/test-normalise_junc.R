context("Testing normalistion of junction counts")

load("example_juncs_w_annot.rda")

juncs <-
  normalise_junc(junc_metadata = example_juncs_w_annot[["metadata"]],
                 raw_count = example_juncs_w_annot[["raw_count"]])

test_that("normalised counts have correct dimensions and class", {
  expect_true(any(class(juncs[["norm_count"]]) == "data.frame"))
  expect_equal(colnames(juncs[["norm_count"]]),
               colnames(juncs[["raw_count"]]))
  expect_equal(nrow(juncs[["norm_count"]]),
               nrow(juncs[["raw_count"]]))
})

test_that("normalised counts have expected values", {
  expect_false(any(is.na(juncs[["norm_count"]])))
  expect_true(all(juncs[["norm_count"]] <= 1))
  expect_true(all(juncs[["norm_count"]] >= 0))
  expect_identical(which(juncs[["raw_count"]] == 0),
                   which(juncs[["norm_count"]] == 0))

  # all clusters with only a single junction should have a value of 1 or 0
  mono_clusters <- which(lengths(mcols(juncs[["metadata"]])[["clusters"]]) == 1)
  expect_true(all(unlist(juncs[["norm_count"]][mono_clusters,]) %in% c(0, 1)))

})

test_that("exon annotation has been correctly retrieved", {

  # test a sample of clusters with >1 junction in them
  # that the clusters have been obtained correctly
  # and the normalisation calculation has been performed correctly
  poly_clusters <- which(lengths(mcols(juncs[["metadata"]])[["clusters"]]) >= 2)

  match <- TRUE

    for(j in sample(poly_clusters, 100)){

      start_hits <- which(start(juncs[["metadata"]][j]) == start(juncs[["metadata"]]))
      end_hits <- which(end(juncs[["metadata"]][j]) == end(juncs[["metadata"]]))
      uniq_hits <- c(start_hits, end_hits) %>%
        unique()

      total_cluster_counts <- lapply(juncs[["raw_count"]][uniq_hits,], sum) %>%
        unlist()

      exp_norm_counts <- juncs[["raw_count"]][j,]/total_cluster_counts
      exp_norm_counts[is.na(exp_norm_counts)] <- 0

      match <- all(match,
                   identical(exp_norm_counts, as.data.frame(juncs[["norm_count"]][j,])))

    }

  expect_true(match)

})

