#' Wrapper for processing coverage
#'
#' `coverage_process` wraps all "coverage_" prefixed functions in
#' [dasper]. This is designed to simplify processing of the coverage data
#' for those familiar or uninterested with the intermediates.
#'
#' @inheritParams coverage_norm
#' @inheritParams coverage_score
#'
#' @return
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   object containing junction data with coverage scores.
#'
#' @examples
#'
#' if (!exists("ref")) {
#'     ref <- "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
#'     ref <- GenomicFeatures::makeTxDbFromGFF(ref)
#' }
#'
#' if (!exists("junctions_processed")) {
#'     junctions_processed <-
#'         junction_process(
#'             junctions_example,
#'             ref,
#'             count_thresh = c("raw" = 5),
#'             n_samp = c("raw" = 1),
#'             width_range = c(25, 1000000),
#'             types = c("ambig_gene", "unannotated"),
#'         )
#' }
#'
#' url <- recount::download_study(
#'     project = "SRP012682",
#'     type = "samples",
#'     download = FALSE
#' )
#' bw_path <- file.path(tempdir(), basename(url[1]))
#'
#' if (!file.exists(bw_path)) {
#'     download.file(url[1], bw_path)
#' }
#'
#' if (!exists("junctions_w_coverage")) {
#'     junctions_w_coverage <-
#'         coverage_process(
#'             junctions_processed,
#'             ref,
#'             unannot_width = 20,
#'             coverage_paths_case = rep(bw_path, 2),
#'             coverage_paths_control = rep(bw_path, 3),
#'             norm_const = 2,
#'             score_func = .zscore,
#'             sd_const = 0.02
#'         )
#' }
#' junctions_w_coverage
#' @family coverage
#' @export
coverage_process <- function(
    junctions,
    ref,
    unannot_width = 20,
    coverage_paths_case,
    coverage_paths_control,
    coverage_chr_control = NULL,
    load_func = .coverage_load,
    bp_param = BiocParallel::SerialParam(),
    norm_const = 1,
    score_func = .zscore,
    ...) {
    print("# Loading and normalising coverage ---------------------------------------------")

    coverage <- coverage_norm(junctions,
        ref,
        unannot_width = 20,
        coverage_paths_case = coverage_paths_case,
        coverage_paths_control = coverage_paths_control,
        coverage_chr_control = NULL,
        load_func = load_func,
        bp_param = bp_param,
        norm_const = norm_const
    )

    print("# Scoring coverage ---------------------------------------------")

    junctions <- coverage_score(junctions, coverage, score_func, ...)

    return(junctions)
}
