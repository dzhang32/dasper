context("Test utility functions")

##### .chr_filter #####

x <- dplyr::tibble(chr = c("1", "2", "3"))

test_that(".chr_filter has correct output", {
    expect_true(tibble::is_tibble(.chr_filter(x, "1")))
    expect_identical(.chr_filter(x, c("1", "2"))[["chr"]], c("1", "2"))
})

test_that(".chr_filter catches user-input errors", {
    expect_error(.chr_filter(x, "chr1"), "No chromosomes in chrs match those of the junction data.")
    expect_warning(.chr_filter(x, c("3", "4", "MT")), "junction data: 4, MT")
})

##### .get_start_end #####

x <- GRanges(c("chr1:1-100", "chr2:2-200"))

x_start_end <- .get_start_end(x)

test_that(".get_start_end has correct output", {
    expect_identical(start(x_start_end$start), start(x))
    expect_identical(end(x_start_end$start), start(x))
    expect_identical(start(x_start_end$end), end(x))
    expect_identical(end(x_start_end$end), end(x))
})

##### .regroup #####

x <- c("A", "B", "C", "D")
groups <- c("1", "1", "3", "4")
all_groups <- 1:5 %>% as.character()
x_regrouped <- .regroup(x, groups, all_groups)

test_that(".regroup has correct output", {
    expect_match(class(x_regrouped), "list")
    expect_identical(x_regrouped[["1"]], c("A", "B"))
    expect_identical(length(x_regrouped[["2"]]), 0L)
    expect_identical(x_regrouped[["3"]], "C")
    expect_identical(x_regrouped[["4"]], "D")
    expect_identical(length(x_regrouped[["5"]]), 0L)
})

##### .merge_CharacterList #####

x <- CharacterList(c("A", "B"), "1")
y <- CharacterList(c("C"), c("2", "3", "4"))

test_that(".merge_CharacterList has correct output", {
    expect_match(class(.merge_CharacterList(x, y)), "CharacterList")
    expect_equal(.merge_CharacterList(x, y)[[1]], c("A", "B", "C"))
    expect_equal(.merge_CharacterList(x, y)[[2]], c("1", "2", "3", "4"))
})

test_that(".merge_CharacterList catches user-input errors", {
    expect_error(
        .merge_CharacterList(x, CharacterList("1")),
        "lengths of x and y should be identical!"
    )
})
