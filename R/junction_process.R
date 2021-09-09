#' Processing junctions
#'
#' @description The set of functions prefixed with "junction_" are used to
#'   process junction data. They are designed to be run in a sequential manner
#'   in the order `junction_annot`, `junction_filter`, `junction_norm`,
#'   `junction_score`. Or, alternatively the wrapper function `junction_process`
#'   can be used to run all 4 of the functions stated above in one go. For more
#'   details of the individual functions, see "Details".
#'
#' @details `junction_process` wraps all "junction_" prefixed functions in
#'   [dasper][dasper::dasper] except [junction_load]. This is designed to
#'   simplify processing of the junction data for those familiar or uninterested
#'   with the intermediates.
#'
#'   `junction_annot` annotates junctions by 1. whether their start and/or end
#'   position precisely overlaps with an annotated exon boundary and 2. whether
#'   that junction matches an intron definition from existing annotation. Using
#'   this information along with the strand, junctions are categorised into
#'   "annotated", "novel_acceptor", "novel_donor", "novel_combo",
#'   "novel_exon_skip", "ambig_gene" and "unannotated".
#'
#'   `junction_filter` filters out "noisy" junctions based on counts, the width
#'   of junctions, annotation category of the junction returned from
#'   [junction_annot] and whether the junction overlaps with a set of
#'   (blacklist) regions.
#'
#'   `junction_norm` normalises the raw junction counts by 1. building junction
#'   clusters by finding junctions that share an acceptor or donor position and
#'   2. calculating a proportion-spliced-in (PSI) for each junction by dividing
#'   the raw junction count by the total number of counts in it's associated
#'   cluster.
#'
#'   `junction_score` will use the counts contained within the "norm"
#'   [assay][SummarizedExperiment::SummarizedExperiment-class] to calculate a
#'   deviation of each patient junction from the expected distribution of
#'   control junction counts. The function used to calculate this abnormality
#'   score can be user-inputted or left as the default z-score. Junctions will
#'   also be labelled based on whether they are up-regulated (+1) or
#'   down-regulated (-1) with respect to controls junction and this information
#'   is stored in the [assay][SummarizedExperiment::SummarizedExperiment-class]
#'   "direction" for use in [outlier_aggregate].
#'
#' @param junctions junction data as a
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'    object.
#' @param ref either path to gtf/gff3 or object of class
#'   [TxDb-class][GenomicFeatures::TxDb-class] or
#'   [EnsDb-class][ensembldb::EnsDb-class].
#'   [EnsDb-class][ensembldb::EnsDb-class] is required if you intend to annotate
#'   junctions with gene symbols/names.
#' @param ref_cols character vector listing the names of the columns in `ref`
#'   for which to annotate junctions with. Must contain "gene_id", used for
#'   categorising junctions.
#' @param ref_cols_to_merge character vector listing which of the annotation
#'   columns `ref_cols` should be merged into in columns to merge into a single
#'   column per junction. Must contain "gene_id", used for categorising
#'   junctions.
#' @param count_thresh named vector with names matching the names of the
#'   [assays][SummarizedExperiment::SummarizedExperiment-class] in `junctions`.
#'   Values denote the number of counts below which a junction will be filtered
#'   out.
#' @param n_samp named vector with names matching the names of the
#'   [assays][SummarizedExperiment::SummarizedExperiment-class] in `junctions`.
#'   Values denotes number of samples that have to express the junction above
#'   the `count_thresh` in order for that junction to not be filtered.
#' @param width_range numeric vector of length 2. The first element denoting the
#'   lower limit of junction width and the second the upper limit. Junctions
#'   with widths outside this range will be filtered out.
#' @param types any junctions matching these types, derived form
#'   [junction_annot] will be filtered out.
#' @param regions any junctions overlapping this set of regions (in a
#'   [GRanges-class][GenomicRanges::GRanges-class] format) will be filtered out.
#' @param score_func function to score junctions by their abnormality. By
#'   default, will use a z-score but can be switched to a user-defined function.
#'   This function must take as input an `x` and `y` argument, containing case
#'   and control counts respectively. This must return a numeric vector equal to
#'   the length of `x` with elements corresponding to a abnormality of each
#'   junction.
#' @param ... additional arguments passed to `score_func`.
#'
#' @return
#' [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object containing filtered, annotated, normalised junction data with
#' abnormality scores.
#'
#' @seealso ENCODE blacklist regions are recommended to be included as `regions`
#'   for `junction_filter` and can be downloaded from
#'   \url{https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz}.
#'    Further information can be found via the publication
#'   \url{https://www.nature.com/articles/s41598-019-45839-z}.
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
#' ##### junction_annot #####
#'
#' junctions <- junction_annot(junctions_example, ref)
#'
#' ##### junction_filter #####
#'
#' junctions <- junction_filter(
#'     junctions,
#'     types = c("ambig_gene", "unannotated")
#' )
#'
#' ##### junction_norm #####
#'
#' junctions <- junction_norm(junctions)
#'
#' ##### junction_score #####
#'
#' junctions <- junction_score(junctions)
#'
#' ##### junction_process #####
#'
#' junctions_processed <- junction_process(
#'     junctions_example,
#'     ref,
#'     types = c("ambig_gene", "unannotated")
#' )
#'
#' # the two objects are equivalent
#' all.equal(junctions_processed, junctions, check.attributes = FALSE)
#' @export
junction_process <- function(junctions,
    ref,
    ref_cols = c("gene_id", "tx_name", "exon_name"),
    ref_cols_to_merge = c("gene_id"),
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

    junctions <- junction_annot(junctions, ref, ref_cols, ref_cols_to_merge)

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
