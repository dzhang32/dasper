context("Test the loading of junctions data")

example_junctions_1_path <- system.file("extdata", "example_junctions_1.txt", package = "dasper", mustWork = TRUE)
example_junctions_2_path <- system.file("extdata", "example_junctions_2.txt", package = "dasper", mustWork = TRUE)

##### .load_STAR #####

example_junctions_1 <- .load_STAR(example_junctions_1_path)

test_that(".load_STAR has the correct output", {
    expect_equal(nrow(example_junctions_1), 10000)
    expect_equal(ncol(example_junctions_1), 5)
    expect_true(tibble::is_tibble(example_junctions_1))
    expect_identical(colnames(example_junctions_1), c("chr", "start", "end", "strand", "count"))
})

##### .junction_merge #####

junctions_all <- .junction_merge(
    junctions_all = NULL,
    junctions = example_junctions_1
)

junctions_all_2 <- .junction_merge(
    junctions_all = junctions_all,
    junctions = example_junctions_1 %>% dplyr::mutate(strand = "*")
)

test_that(".junction_merge has the correct output", {
    expect_true(any(colnames(junctions_all) == "count_1"))
    expect_true(all(c("count_1", "count_2") %in% colnames(junctions_all_2)))
    expect_identical(junctions_all["strand"], junctions_all_2["strand"])
})

##### .control_coord_convert #####

ensembl <- example_junctions_1 %>%
    dplyr::mutate(chr = "MT")
ucsc <- ensembl %>%
    dplyr::mutate(
        chr = stringr::str_c("chr", chr),
        start = start - 1,
        end = end - 1
    )

# control data is currently always UCSC based
ucsc_2 <- .control_coord_convert(ucsc, coord_system = "ucsc")
ensembl_2 <- .control_coord_convert(ucsc, coord_system = "ensembl")

test_that(".control_coord_convert has the correct output", {
    expect_identical(ucsc, ucsc_2)
    expect_equivalent(ensembl, ensembl_2)
})

test_that(".control_coord_convert catches user-input errors", {
    expect_error(
        .control_coord_convert(ucsc, coord_system = "not_a_coord_system"),
        "coord_system must be one of 'ensembl' or 'ucsc'"
    )
})

##### junction_load #####

junctions <-
    junction_load(
        junction_paths = c(example_junctions_1_path)
    )

junctions_w_case <-
    junction_load(
        junction_paths = c(example_junctions_1_path, example_junctions_2_path),
        metadata = dplyr::tibble(samp_id = c("example_1", "example_2")),
        controls = c(TRUE, FALSE)
    )

junctions_w_control <-
    junction_load(
        junction_paths = c(example_junctions_1_path, example_junctions_2_path),
        controls = "fibroblasts",
        chrs = c("21"),
        coord_system = "ensembl"
    )

test_that("junction_load has correct output", {
    expect_match(class(junctions), "RangedSummarizedExperiment")
    expect_identical(SummarizedExperiment::colData(junctions)[["samp_id"]], c("samp_1"))
    expect_identical(SummarizedExperiment::colData(junctions)[["case_control"]], c("case"))
    expect_identical(SummarizedExperiment::assay(junctions) %>% colnames(), c("count_1"))

    expect_match(class(junctions_w_case), "RangedSummarizedExperiment")
    expect_identical(SummarizedExperiment::colData(junctions_w_case)[["case_control"]], c("control", "case"))
    expect_identical(SummarizedExperiment::assay(junctions_w_case) %>% colnames(), c("count_1", "count_2"))
    expect_identical(SummarizedExperiment::colData(junctions_w_case)[["samp_id"]], c("example_1", "example_2"))

    expect_match(class(junctions_w_control), "RangedSummarizedExperiment")
    expect_identical(
        SummarizedExperiment::colData(junctions_w_control)[["case_control"]],
        c(rep("case", 2), rep("control", ncol(SummarizedExperiment::assay(junctions_w_control)) - 2))
    )
    expect_identical(
        SummarizedExperiment::assay(junctions_w_control) %>%
            colnames() %>%
            stringr::str_replace("_.*", "") %>%
            unique(),
        c("count", "gtex")
    )
    expect_identical(
        SummarizedExperiment::rowRanges(junctions_w_control) %>%
            GenomicRanges::seqnames() %>%
            unique() %>%
            as.character(),
        "21"
    )
})

test_that("junction_load catches user-input errors", {
    expect_error(
        junction_load(
            junction_paths = c(example_junctions_1_path, example_junctions_2_path),
            metadata = tibble(samp_id = c("example_1", "example_2")),
            controls = "must_be_lgl"
        ),
        "Controls argument must be a logical vector of same length as the number of junction paths"
    )
    expect_error(
        junction_load(
            junction_paths = c(example_junctions_1_path, example_junctions_2_path),
            metadata = tibble(samp_id = c("example_1", "example_2")),
            controls = TRUE
        ),
        "Controls argument must be a logical vector of same length as the number of junction paths"
    )
    expect_error(
        junction_load(
            junction_paths = c(example_junctions_1_path, example_junctions_2_path),
            metadata = tibble(samp_id = c("example_1", "example_2")),
            controls = "not_a_ctrl"
        ),
        "Controls argument must be a logical vector of same length as the number of junction paths"
    )
})
