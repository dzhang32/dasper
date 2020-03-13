context("Testing common utility functions")

gr <- GRanges((stringr::str_c("chr1:", 1:10, "-", 11:20)))

test_that(".get_gr_for_start_end has correct output", {
  expect_match(class(.get_gr_for_start_end(gr)), "list")
  expect_equal(start(.get_gr_for_start_end(gr)[["start"]]),
               start(gr))
  expect_equal(end(.get_gr_for_start_end(gr)[["start"]]),
               start(gr))
  expect_equal(start(.get_gr_for_start_end(gr)[["end"]]),
               end(gr))
  expect_equal(end(.get_gr_for_start_end(gr)[["end"]]),
               end(gr))
})

x <- CharacterList(c("A", "B"), "1")
y <- CharacterList(c("C"), c("2", "3", "4"))
names(x) <- c("x", "y")
names(y) <- c("x", "y")

test_that("merging lists gives correct output", {
  expect_match(class(.merge_lists(x, y)), "CharacterList")
  expect_equal(.merge_lists(x, y)[["x"]], c("A", "B", "C"))
  expect_equal(.merge_lists(x, y)[["y"]], c("1", "2", "3", "4"))
  expect_error(.merge_lists(x, CharacterList(c("C"), c("2", "3", "4"))),
               "names of x and y lists should be identical!")
})
