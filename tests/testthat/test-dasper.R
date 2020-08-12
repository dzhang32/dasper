context("Testing dasper")

##### reticulate set-up #####

if (.Platform$OS.type != "windows") {

    # force reticulate to use the python3 install
    # if windows skip this step
    reticulate::use_python(Sys.which("python3"), required = TRUE)
}

##### Set up random score data #####

junctions <- junctions_example[, colData(junctions_example)[["case_control"]] == "case"]

# add random scores and direction to save time

direction <- matrix(
    data = sample(c(1, -1), dim(junctions)[1] * dim(junctions)[2], replace = TRUE),
    nrow = dim(junctions)[1],
    ncol = dim(junctions)[2]
)

score <- matrix(
    data = sample(-100:100, dim(junctions)[1] * dim(junctions)[2], replace = TRUE),
    nrow = dim(junctions)[1],
    ncol = dim(junctions)[2]
)

coverage_score <- matrix(
    data = sample(-100:100, dim(junctions)[1] * dim(junctions)[2], replace = TRUE),
    nrow = dim(junctions)[1],
    ncol = dim(junctions)[2]
)

colnames(direction) <- dimnames(junctions)[[2]]
colnames(score) <- dimnames(junctions)[[2]]
colnames(coverage_score) <- dimnames(junctions)[[2]]

assays(junctions)[["direction"]] <- direction
assays(junctions)[["score"]] <- score
assays(junctions)[["coverage_score"]] <- coverage_score

##### outlier_detect #####

junctions_w_outlier_score <- outlier_detect(junctions,
    feature_names = c("score", "coverage_score"),
    random_state = 32L
)

up_indexes <- which(assays(junctions)[["direction"]][, 1] == 1)
outlier_up <- .outlier_score(
    features =
        data.frame(
            score = assays(junctions)[["score"]][, 1][up_indexes],
            coverage_score = assays(junctions)[["coverage_score"]][, 1][up_indexes]
        ),
    random_state = 32L
)

down_indexes <- which(assays(junctions)[["direction"]][, 1] == -1)
outlier_down <- .outlier_score(
    features =
        data.frame(
            score = assays(junctions)[["score"]][, 1][down_indexes],
            coverage_score = assays(junctions)[["coverage_score"]][, 1][down_indexes]
        ),
    random_state = 32L
)

test_that("outlier_detect has the correct output", {

    # if junctions have been reordered either one of the direction
    # or scores of up/down would not match
    expect_identical(
        assays(junctions_w_outlier_score)[["direction"]],
        assays(junctions)[["direction"]]
    )

    expect_identical(
        outlier_up %>% as.numeric(),
        assays(junctions_w_outlier_score)[["outlier_score"]][, 1][up_indexes]
    )

    expect_identical(
        outlier_down %>% as.numeric(),
        assays(junctions_w_outlier_score)[["outlier_score"]][, 1][down_indexes]
    )
})

test_that("outlier_detect catches user-input errors", {
    expect_error(
        outlier_detect(junctions, feature_names = c("not_an_assay")),
        "Assays does not contain the following: "
    )

    assays(junctions)[["direction"]] <- NULL

    expect_error(
        outlier_detect(junctions),
        "junctions must contain a 'direction' assay"
    )
})
