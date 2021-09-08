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

##### .coverage_load #####

megadepth::install_megadepth()

# obtain path to example bw on recount2
url <- recount::download_study(
    project = "SRP012682",
    type = "samples",
    download = FALSE
)

bw_path <- dasper:::.file_cache(url[1])

# take 5 junctions from the end and 5 from the top to test in that order
# to make sure order returned by .coverage_load() matches input ranges rather than
# chromosome order
junctions_to_use <- c(length(junctions_example):(length(junctions_example) - 4), 1:5)
junctions <- junctions_example[junctions_to_use] %>%
    rowRanges()

junctions_sorted <- junctions %>% sort()

mcols(junctions)[["coverage"]] <-
    .coverage_load(
        regions = junctions,
        coverage_path = bw_path,
        sum_fun = "mean",
        chr_format = "chr"
    )

mcols(junctions_sorted)[["coverage"]] <-
    .coverage_load(
        regions = junctions_sorted,
        coverage_path = bw_path,
        sum_fun = "mean",
        chr_format = "chr"
    )

coverage_rt <-
    .coverage_load(
        regions = junctions_sorted,
        coverage_path = bw_path,
        sum_fun = "mean",
        chr_format = "chr",
        method = "rt"
    ) %>%
    lapply(FUN = mean) %>%
    unlist() %>%
    round(2)

test_that(".coverage_load has correct output", {

    # make sure the order of returned coverage is same as inputted regions
    expect_identical(
        mcols(sort(junctions))[["coverage"]],
        mcols(junctions_sorted)[["coverage"]]
    )

    expect_equal(
        mcols(junctions_sorted)[["coverage"]],
        coverage_rt
    )
})

##### .chr_check #####

x <- GRanges(c("chr1:1-100", "chr2:2-200"))
y <- GRanges(c("1:1-100", "2:2-200"))

test_that(".chr_check has correct output", {
    expect_identical(x, .chr_check(x, "chr"))
    expect_identical(y, .chr_check(x, "no_chr"))
    expect_identical(x, .chr_check(y, "chr"))
    expect_identical(y, .chr_check(y, "no_chr"))
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

##### .outlier_score #####

features <- data.frame(
    index = 1:10,
    feat_1 = 1:10,
    feat_2 = 1:10
)

features_desc <- features %>%
    dplyr::arrange(desc(index))

features[["score"]] <- .outlier_score(features, random_state = 32L)
features_desc[["score"]] <- .outlier_score(features_desc, random_state = 32L)

test_that(".outlier_score has correct output", {

    # check order not changed
    expect_identical(
        features[["score"]],
        features_desc %>% dplyr::arrange(index) %>% .[["score"]]
    )
})

##### ref_load #####

ref <- GenomicState::GenomicStateHub(version = "31", genome = "hg38", filetype = "TxDb")[[1]]

ref <- ref_load(ref)

test_that("ref_load has correct output", {
    expect_equal(GenomicFeatures::organism(ref), "Homo sapiens")
})

test_that("ref_load catches user input errors", {
    expect_error(
        ref_load(2),
        "ref must either be a path to the .gtf/gff3 file "
    )
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
