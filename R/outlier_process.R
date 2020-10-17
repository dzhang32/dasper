#' Processing outliers
#'
#' @description The set of functions prefixed with "outlier_" are used to detect
#'   outliers. They are designed to be run after you have extracted your
#'   junctions and coverage based features, in the order `outlier_detect`,
#'   `outlier_aggregate`. Or, alternatively the wrapper function
#'   `outlier_process` can be used to run the 2 functions stated above in one
#'   go. For more details of the individual functions, see "Details".
#'
#' @details `outlier_process` wraps all "outlier_" prefixed functions in
#'   [dasper]. This is designed to simplify processing of the detecting outlier
#'   junctions for those familiar or uninterested with the intermediates.
#'
#'   `outlier_detect` will use the features in
#'   [assays][SummarizedExperiment::SummarizedExperiment-class] named
#'   `feature_names` as input into an unsupervised outlier detection algorithm
#'   to score each junction based on how outlier-y it looks in relation to other
#'   junctions in the patient. The default expected `score` and `coverage_score`
#'   features can be calculated using the [junction_process] and
#'   [coverage_process] respectively.
#'
#'   `outlier_aggregate` will aggregate the outlier scores into a cluster-level.
#'   It will then rank each cluster based on this aggregated score and annotate
#'   each cluster with it's associated gene and transcript.
#'
#' @inheritParams junction_annot
#'
#' @param feature_names names of assays in `junctions` that are to be used as
#'   input into the outlier detection model.
#' @param bp_param a
#'   [BiocParallelParam-class][BiocParallel::BiocParallelParam-class] instance
#'   denoting whether to parallelise the calculating of outlier scores across
#'   samples.
#' @param ... additional arguments passed to the outlier detection model
#'   (isolation forest) for setting parameters.
#' @param samp_id_col name of the column in the
#'   [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class]
#'   that details the sample ids.
#'
#' @return `DataFrame` with one row per cluster detailing each cluster's
#'   associated junctions, outlier scores, ranks and genes.
#'
#' @seealso for more details on the isolation forest model used:
#'   https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.IsolationForest.html
#'
#' @examples
#'
#' ##### Set up txdb #####
#'
#' # use GenomicState to load txdb (GENCODE v31)
#' ref <- GenomicState::GenomicStateHub(
#'     version = "31",
#'     genome = "hg38",
#'     filetype = "TxDb"
#' )[[1]]
#'
#' ##### Set up BigWig #####
#'
#' # obtain path to example bw on recount2
#' bw_path <- recount::download_study(
#'     project = "SRP012682",
#'     type = "samples",
#'     download = FALSE
#' )[[1]]
#' \dontshow{
#' # cache the bw for speed in later
#' # examples/testing during R CMD Check
#' bw_path <- dasper:::.file_cache(bw_path)
#' }
#'
#' ##### junction_process #####
#'
#' junctions_processed <- junction_process(
#'     junctions_example,
#'     ref,
#'     types = c("ambig_gene", "unannotated"),
#' )
#'
#' ##### coverage_process #####
#'
#' junctions_w_coverage <- coverage_process(
#'     junctions_processed,
#'     ref,
#'     coverage_paths_case = rep(bw_path, 2),
#'     coverage_paths_control = rep(bw_path, 3)
#' )
#' \donttest{
#' ##### outlier_process #####
#'
#' # this wrapper will obtain outlier scores identical to those
#' # obtained through running the individual wrapped functions shown below
#' outlier_process(junctions_w_coverage)
#' }
#'
#' ##### outlier_detect #####
#'
#' junctions_w_outliers <- outlier_detect(junctions_w_coverage)
#'
#' ##### outlier_aggregate #####
#'
#' outlier_scores <- outlier_aggregate(junctions_w_outliers)
#' @export
outlier_process <- function(
    junctions,
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
