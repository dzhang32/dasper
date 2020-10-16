#' Wrapper for processing junctions
#'
#' `junction_process` wraps all "junction_" prefixed functions in
#' [dasper][dasper::dasper] except [junction_load]. This is designed to simplify
#' processing of the junction data for those familiar or uninterested with the
#' intermediates.
#'
#' @inheritParams junction_filter
#' @inheritParams junction_annot
#' @inheritParams junction_norm
#' @inheritParams junction_score
#'
#' @return
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   object containing filtered, annotated, normalised junction data with abnormality scores.
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
#' @family junction
#' @export
junction_process <- function(
    junctions,
    ref,
    count_thresh = c("raw" = 5),
    n_samp = c("raw" = 1),
    width_range = NULL,
    types = NULL,
    regions = NULL,
    score_func = .zscore,
    ...) {
    if (!is.null(count_thresh) | !is.null(width_range) | is.null(regions)) {
        print("# Filtering junctions -----------------------------------------------------")

        # first filter by count/width/overlap to save time for junction_annot
        junctions <- junction_filter(junctions,
            count_thresh = count_thresh,
            n_samp = n_samp,
            width_range = width_range,
            regions = regions
        )
    }

    print("# Annotating junctions ----------------------------------------------------")

    junctions <- junction_annot(junctions, ref)

    if (!is.null(types)) {
        print("# Filtering junctions -----------------------------------------------------")

        junctions <- junction_filter(junctions,
            count_thresh = NULL,
            n_samp = NULL,
            types = types
        )
    }

    print("# Normalise junctions -----------------------------------------------------")

    junctions <- junction_norm(junctions)

    print("# Score junctions ---------------------------------------------------------")

    junctions <- junction_score(junctions, score_func, ...)

    return(junctions)
}
