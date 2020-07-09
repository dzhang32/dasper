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

junction_df_all <- .junction_merge(
    junction_df_all = NULL,
    example_junctions_1
)

junction_df_all_2 <- .junction_merge(
    junction_df_all,
    example_junctions_1 %>% dplyr::mutate(strand = "*")
)

test_that(".junction_merge has the correct output", {
    expect_true(any(colnames(junction_df_all) == "count_1"))
    expect_true(all(c("count_1", "count_2") %in% colnames(junction_df_all_2)))
    expect_identical(junction_df_all["strand"], junction_df_all_2["strand"])
})

##### junction_load #####

junctions_control <-
    junction_load(
        junction_paths = c(example_junctions_1_path),
        metadata = dplyr::tibble(samp_id = c("example_1"))
    )

junctions_w_case <-
    junction_load(
        junction_paths = c(example_junctions_1_path, example_junctions_2_path),
        metadata = dplyr::tibble(samp_id = c("example_1", "example_2")),
        controls = c(TRUE, FALSE)
    )

test_that("junction_load has correct output", {
    expect_match(class(junctions_control), "RangedSummarizedExperiment")
    expect_identical(SummarizedExperiment::colData(junctions_control)[["case_control"]], c("case"))
    expect_identical(SummarizedExperiment::colData(junctions_w_case)[["case_control"]], c("control", "case"))
    expect_identical(SummarizedExperiment::assay(junctions_control) %>% colnames(), c("count_1"))
    expect_identical(SummarizedExperiment::assay(junctions_w_case) %>% colnames(), c("count_1", "count_2"))
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
})
