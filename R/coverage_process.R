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
#' # use Genomic state to load txdb (GENCODE v31)
#' ref <- GenomicState::GenomicStateHub(version = "31", genome = "hg38", filetype = "TxDb")[[1]]
#'
#' junctions_processed <- junction_process(
#'     junctions_example,
#'     ref,
#'     count_thresh = c("raw" = 5),
#'     n_samp = c("raw" = 1),
#'     types = c("ambig_gene", "unannotated"),
#' )
#'
#' # obtain path to example bw on recount2
#' url <- recount::download_study(
#'     project = "SRP012682",
#'     type = "samples",
#'     download = FALSE
#' )
#'
#' bw_path <- dasper:::.file_cache(url[1])
#'
#' junctions_w_coverage <- coverage_process(
#'     junctions_processed,
#'     ref,
#'     coverage_paths_case = rep(bw_path, 2),
#'     coverage_paths_control = rep(bw_path, 3)
#' )
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
