context("Loading and merging junctions")

example_juncs_1 <- .load_STAR("../../data-raw/example_juncs_1.txt", sample_id = "example_1")

juncs <-
  merge_junc(junc_paths = c("../../data-raw/example_juncs_1.txt",
                            "../../data-raw/example_juncs_2.txt"),
             sample_ids = c("eg1", "eg2"),
             load_func = .load_STAR,
             chr_to_filter = NULL)

test_that(".load_STAR has the correct output", {
  expect_equal(nrow(example_juncs_1), 10000)
  expect_equal(ncol(example_juncs_1), 5)
  expect_true(tibble::is_tibble(example_juncs_1))
})

test_that("merging juncs has correct output", {
  expect_gte(nrow(juncs$raw_count), 10000)
  expect_true(is.list(juncs))
  expect_match(class(juncs$metadata), "GRanges")
  expect_true(tibble::is_tibble(juncs$raw_count))
  expect_false(any(is.na(juncs)))
})

test_that("merging juncs has catches errors", {

  expect_error(merge_junc(junc_paths = c("one", "two", "three"),
                          sample_ids = c("one", "two")),
               "Number of junc_paths does not equal to the number of sample_ids")

  expect_error(merge_junc(junc_paths = c("../../data-raw/example_juncs_1.txt",
                                         "../../data-raw/example_juncs_2.txt"),
                          sample_ids = c("eg1", "eg2"),
                          chr_to_filter = "not_a_chr"),
               "No chromosomes in chr_to_filter match those of junction data")

})
