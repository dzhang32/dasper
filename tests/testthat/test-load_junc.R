context("Test the loading of junctions data")

example_juncs_1_path <- system.file("extdata", "example_juncs_1.txt", package = "dasper", mustWork = TRUE)
example_juncs_2_path <- system.file("extdata", "example_juncs_2.txt", package = "dasper", mustWork = TRUE)

##### .load_STAR #####

example_juncs_1 <- .load_STAR(example_juncs_1_path)

test_that(".load_STAR has the correct output", {
    expect_equal(nrow(example_juncs_1), 10000)
    expect_equal(ncol(example_juncs_1), 5)
    expect_true(tibble::is_tibble(example_juncs_1))
    expect_identical(colnames(example_juncs_1), c("chr", "start", "end", "strand", "count"))
})

##### .junc_merge #####

junc_df_all <- .junc_merge(
    junc_df_all = tibble(),
    example_juncs_1,
    i = 1
)

junc_df_all_2 <- .junc_merge(junc_df_all,
    example_juncs_1 %>% dplyr::mutate(strand = "*"),
    i = 2
)

test_that(".junc_merge has the correct output", {
    expect_true(any(colnames(junc_df_all) == "count_1"))
    expect_true(all(c("count_1", "count_2") %in% colnames(junc_df_all_2)))
    expect_identical(junc_df_all["strand"], junc_df_all_2["strand"])
})

##### junc_load #####

juncs_control <-
    junc_load(
        junc_paths = c(example_juncs_1_path),
        metadata = dplyr::tibble(samp_id = c("example_1"))
    )

juncs_w_case <-
    junc_load(
        junc_paths = c(example_juncs_1_path, example_juncs_2_path),
        metadata = dplyr::tibble(samp_id = c("example_1", "example_2")),
        controls = c(TRUE, FALSE)
    )

test_that("junc_load has correct output", {
    expect_match(class(juncs_control), "RangedSummarizedExperiment")
    expect_identical(SummarizedExperiment::colData(juncs_control)[["case_control"]], c("case"))
    expect_identical(SummarizedExperiment::colData(juncs_w_case)[["case_control"]], c("control", "case"))
})

test_that("junc_load catches user-input errors", {
    expect_error(
        junc_load(
            junc_paths = c(example_juncs_1_path, example_juncs_2_path),
            metadata = tibble(samp_id = c("example_1", "example_2")),
            controls = "must_be_lgl"
        ),
        "Controls argument must be a logical vector of same length as the number of junction paths"
    )
    expect_error(
        junc_load(
            junc_paths = c(example_juncs_1_path, example_juncs_2_path),
            metadata = tibble(samp_id = c("example_1", "example_2")),
            controls = TRUE
        ),
        "Controls argument must be a logical vector of same length as the number of junction paths"
    )
})
