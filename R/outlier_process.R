#' Wrapper for processing outliers
#'
#' `outlier_process` wraps all "outlier_" prefixed functions in
#' [dasper]. This is designed to simplify processing of the detecting
#' outlier junctions for those familiar or uninterested with the intermediates.
#'
#' @inheritParams junction_annot
#' @inheritParams outlier_detect
#' @inheritParams outlier_aggregate
#'
#' @return `DataFrame` with one row per cluster detailing each cluster's
#'   associated junctions, outlier scores, ranks and genes.
#'
#' @family outlier
#' @export
#'
#' @examples
#'
#' if (.Platform$OS.type != "windows") {
#'     # tell reticulate to use the python3 install
#'     # if windows skip this step
#'     reticulate::use_python(Sys.which("python3"), required = TRUE)
#' }
#'
#' if (!exists("ref")) {
#'     # use Genomic state to load txdb (GENCODE v31)
#'     ref <- GenomicState::GenomicStateHub(version = "31", genome = "hg38", filetype = "TxDb")[[1]]
#'     # convert seqlevels to match junctions
#'     seqlevels(ref) <- stringr::str_replace(seqlevels(ref), "chr", "")
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
#'     GenomeInfoDb::seqlevels(junctions_processed) <-
#'         paste0("chr", GenomeInfoDb::seqlevels(junctions_processed))
#' }
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
#' if (!exists("junctions_w_coverage")) {
#'     junctions_w_coverage <-
#'         coverage_process(
#'             junctions_processed,
#'             ref,
#'             coverage_paths_case = rep(bw_path, 2),
#'             coverage_paths_control = rep(bw_path, 3)
#'         )
#' }
#'
#' outlier_process(junctions_w_coverage)
outlier_process <- function(junctions,
    feature_names = c("score", "coverage_score"),
    samp_id_col = "samp_id",
    bp_param = BiocParallel::SerialParam(),
    ...) {
    print("# Detecting outliers using an isolation forest ------------------------------")

    junctions <- outlier_detect(
        junctions = junctions,
        feature_names = feature_names,
        bp_param = bp_param,
        ...
    )

    print("# Aggregating outlier scores to cluster-level -----------------------------")

    outlier_scores_tidy <- outlier_aggregate(
        junctions = junctions,
        samp_id_col = samp_id_col,
        bp_param = bp_param
    )

    return(outlier_scores_tidy)
}
