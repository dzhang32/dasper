context("Testing junction filtering")

load("example_juncs_w_annot.rda")

juncs <-
  filter_junc(junc_metadata = example_juncs_w_annot[["metadata"]],
              raw_count = example_juncs_w_annot[["raw_count"]],
              count_thresh = 5,
              n_samp = 1,
              junc_width = c(20, 500000),
              junc_cats = c("ambig_gene", "none"),
              blacklist_gr = GRanges("1:1-100000000"))

test_that("incorrect junction categories are picked up", {
  expect_warning(filter_junc(junc_metadata = example_juncs_w_annot[["metadata"]],
                             raw_count = example_juncs_w_annot[["raw_count"]],
                             junc_cats = c("ambig_gene", "none", "check", "check2")),
                 regexp = "The following junction categories are not expected so ignored")
})

test_that("each filter has been applied correctly", {

  # first check that without filters applied
  # there are junctions that match the filter criteria
  count_test <-
  example_juncs_w_annot[["raw_count"]] %>%
    apply(MARGIN = 1, FUN = function(x) (sum(x >= 5)))

  expect_true(any(count_test < 1))
  expect_true(any(width(example_juncs_w_annot[["metadata"]]) < 20))
  expect_true(any(width(example_juncs_w_annot[["metadata"]]) > 500000))
  expect_true(any(example_juncs_w_annot[["metadata"]]$junc_cat %in% c("ambig_gene", "none")))
  expect_true(length(findOverlaps(example_juncs_w_annot$metadata, GRanges("1:1-100000000"))) > 0)

  # then check none are left after filtering
  count_test <-
  juncs[["raw_count"]] %>%
    apply(MARGIN = 1, FUN = function(x) (sum(x >= 5)))

  expect_true(all(count_test >= 1))
  expect_false(any(width(juncs$metadata) < 20))
  expect_false(any(width(juncs$metadata) > 1e6))
  expect_false(any(juncs$metadata$junc_cat %in% c("ambig_gene", "none")))
  expect_true(length(findOverlaps(juncs$metadata, GRanges("1:1-100000000"))) == 0)

})
