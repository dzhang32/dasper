context("Testing common helper functions")

gr <- GenomicRanges::GRanges((stringr::str_c("chr1:", 1:10, "-", 11:20)))

test_that(".get_gr_for_start_end has correct output", {

  expect_match(class(.get_gr_for_start_end(gr)), "list")

  expect_equal(GenomicRanges::start(.get_gr_for_start_end(gr)[["start"]]),
               GenomicRanges::start(gr))
  expect_equal(GenomicRanges::end(.get_gr_for_start_end(gr)[["start"]]),
               GenomicRanges::start(gr))

  expect_equal(GenomicRanges::start(.get_gr_for_start_end(gr)[["end"]]),
               GenomicRanges::end(gr))
  expect_equal(GenomicRanges::end(.get_gr_for_start_end(gr)[["end"]]),
               GenomicRanges::end(gr))

})
