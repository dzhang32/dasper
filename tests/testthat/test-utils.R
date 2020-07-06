context("Test common utility functions")

##### .chr_filter #####

x <- dplyr::tibble(chr = c("1", "2", "3"))
.chr_filter(x, "1")

test_that(".chr_filter has the correct output", {
    expect_true(tibble::is_tibble(.chr_filter(x, "1")))
    expect_identical(.chr_filter(x, c("1", "2"))[["chr"]], c("1", "2"))
})

test_that(".chr_filter catches user-input errors", {
    expect_error(.chr_filter(x, "chr1"), "No chromosomes in chrs match those of the junction data.")
    expect_warning(.chr_filter(x, c("3", "4", "MT")), "junction data: 4, MT")
})
