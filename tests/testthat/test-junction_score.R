context("Testing junction scoring")

##### junction_filter #####

junctions <- junction_norm(junctions_example)
junctions_w_score <- junction_score(junctions,
    score_func = .zscore,
    sd_const = 0.02
) # try an sd_const other than default

test_that("junction_score has correct output", {
    expect_true(is(junctions_w_score, "RangedSummarizedExperiment"))
    expect_identical(dim(junctions_w_score)[1], dim(junctions)[1])
    expect_identical(dim(junctions_w_score)[2], sum(colData(junctions)[["case_control"]] == "case"))
    expect_identical(names(assays(junctions_w_score)), c("raw", "norm", "score"))
})

junctions_no_colData <- junctions_example
SummarizedExperiment::colData(junctions_no_colData)[["case_control"]] <- NULL

test_that("junction_score catches user-input errors", {
    expect_error(
        junction_score(junctions_no_colData),
        "must include the column 'case_control'"
    )

    expect_error(
        junction_score(junctions_example),
        "Junctions must include the 'norm' assay"
    )
})

# vectorised method of performing z-score
# test whether this produces equivalent output to junction_score
case_count <- assays(junctions)[["norm"]][, colData(junctions)[["case_control"]] == "case"]
control_count <- assays(junctions)[["norm"]][, colData(junctions)[["case_control"]] == "control"]

control_mean <-
    control_count %>%
    apply(MARGIN = 1, FUN = mean)

control_sd <-
    control_count %>%
    apply(MARGIN = 1, FUN = sd)

case_score <- (case_count - control_mean) / (control_sd + 0.02)

test_that("zscore has been calculated correctly", {
    expect_equivalent(
        assays(junctions_w_score)[["score"]] %>% unlist(),
        case_score %>% unlist()
    )
})
