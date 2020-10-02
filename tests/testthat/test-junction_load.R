context("Test the loading of junctions data")

junctions_example_1_path <- system.file("extdata", "junctions_example_1.txt", package = "dasper", mustWork = TRUE)
junctions_example_2_path <- system.file("extdata", "junctions_example_2.txt", package = "dasper", mustWork = TRUE)

##### .STAR_load #####

junctions_example_1 <- .STAR_load(junctions_example_1_path)

test_that(".STAR_load has the correct output", {
    expect_equal(nrow(junctions_example_1), 10000)
    expect_equal(ncol(junctions_example_1), 5)
    expect_true(is(junctions_example_1, "data.frame"))
    expect_identical(colnames(junctions_example_1), c("chr", "start", "end", "strand", "count"))
})

##### .junction_merge #####

junctions_all <- .junction_merge(
    junctions_all = NULL,
    junctions = junctions_example_1
)

junctions_all_2 <- .junction_merge(
    junctions_all = junctions_all,
    junctions = junctions_example_1 %>% dplyr::mutate(strand = "*")
)

test_that(".junction_merge has the correct output", {
    expect_true(any(colnames(junctions_all) == "count_1"))
    expect_true(all(c("count_1", "count_2") %in% colnames(junctions_all_2)))
    expect_identical(junctions_all["strand"], junctions_all_2["strand"])
})

##### .control_coord_convert #####

# control data (raw GTEx data downloaded via snaptron)
# uses 1-based co-ordinates, yet prefixes chromsomes with a "chr"
control <- junctions_example_1 %>%
    dplyr::mutate(chr = chr %>% stringr::str_c("chr", .) %>%
        stringr::str_replace("MT", "M"))

# generate ensembl/ucsc example junctions
ensembl <- junctions_example_1
ucsc <- junctions_example_1 %>%
    dplyr::mutate(
        chr = stringr::str_c("chr", chr) %>%
            stringr::str_replace("MT", "M"),
        start = start - 1,
        end = end - 1
    )

ucsc_2 <- .control_coord_convert(control, coord_system = "ucsc")
ensembl_2 <- .control_coord_convert(control, coord_system = "ensembl")

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
        junction_paths = c(junctions_example_1_path)
    )

junctions_w_case <-
    junction_load(
        junction_paths = c(junctions_example_1_path, junctions_example_2_path),
        metadata = dplyr::tibble(samp_id = c("example_1", "example_2")),
        controls = c(TRUE, FALSE)
    )

junctions_w_control <-
    junction_load(
        junction_paths = c(junctions_example_1_path, junctions_example_2_path),
        controls = "fibroblasts",
        chrs = c("21", "MT"),
        coord_system = "ensembl"
    )

test_that("junction_load has correct output", {
    expect_match(class(junctions), "RangedSummarizedExperiment")
    expect_identical(colData(junctions)[["samp_id"]], c("samp_1"))
    expect_identical(colData(junctions)[["case_control"]], c("case"))
    expect_identical(assays(junctions)[["raw"]] %>% colnames(), c("count_1"))

    expect_match(class(junctions_w_case), "RangedSummarizedExperiment")
    expect_identical(colData(junctions_w_case)[["case_control"]], c("control", "case"))
    expect_identical(assays(junctions_w_case)[["raw"]] %>% colnames(), c("count_1", "count_2"))
    expect_identical(colData(junctions_w_case)[["samp_id"]], c("example_1", "example_2"))

    expect_match(class(junctions_w_control), "RangedSummarizedExperiment")
    expect_identical(
        colData(junctions_w_control)[["case_control"]],
        c(rep("case", 2), rep("control", ncol(assays(junctions_w_control)[["raw"]]) - 2))
    )
    expect_identical(
        assays(junctions_w_control)[["raw"]] %>%
            colnames() %>%
            stringr::str_replace("_.*", "") %>%
            unique(),
        c("count", "gtex")
    )
    expect_identical(
        rowRanges(junctions_w_control) %>%
            GenomicRanges::seqnames() %>%
            unique() %>%
            as.character(),
        c("21", "MT")
    )
})

test_that("junction_load catches user-input errors", {
    expect_error(
        junction_load(
            junction_paths = c(junctions_example_1_path, junctions_example_2_path),
            metadata = tibble(samp_id = c("example_1", "example_2")),
            controls = "must_be_lgl"
        ),
        "Controls argument must be a logical vector of same length as the number of junction paths"
    )
    expect_error(
        junction_load(
            junction_paths = c(junctions_example_1_path, junctions_example_2_path),
            metadata = tibble(samp_id = c("example_1", "example_2")),
            controls = TRUE
        ),
        "Controls argument must be a logical vector of same length as the number of junction paths"
    )
    expect_error(
        junction_load(
            junction_paths = c(junctions_example_1_path, junctions_example_2_path),
            metadata = tibble(samp_id = c("example_1", "example_2")),
            controls = "not_a_ctrl"
        ),
        "Controls argument must be a logical vector of same length as the number of junction paths"
    )
})
