#' Check chromosomes are correctly formatted
#'
#' `.chr_filter` will check whether the object `x` is matches the
#' format of chr_format. If not, will convert to the chromosome format.
#'
#' @param x object of [GRangesList-class][GenomicRanges::GRangesList-class] or
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   class.
#' @param chr_format one of "chr" or "no_chr".
#'
#' @return `x` with chromosomes matching `chr_format`.
#'
#' @keywords internal
#' @noRd
.chr_check <- function(x, chr_format) {
    if (chr_format == "chr" & !any(stringr::str_detect(GenomeInfoDb::seqlevels(x), "chr"))) {
        GenomeInfoDb::seqlevels(x) <- stringr::str_c("chr", GenomeInfoDb::seqlevels(x))
    } else if (chr_format == "no_chr" & any(stringr::str_detect(seqnames(x), "chr"))) {
        GenomeInfoDb::seqlevels(x) <- GenomeInfoDb::seqlevels(x) %>%
            stringr::str_replace("chr", "")
    }

    return(x)
}

#' Filter a data.frame for chromosomes to keep
#'
#' `.chr_filter` will filter a data.frame using [filter][dplyr::filter], keeping
#' only chromosomes listed in the chrs argument. It includes a couple of checks
#' to notify users whether some or all chrs are not found in their data. This
#' aims to catch situations where alternative nomenclature may be used - e.g.
#' "M" vs "MT", or "chr1" vs "1".
#'
#' @param x data.frame containing at minimum a column named "chr".
#' @inheritParams junction_load
#'
#' @return data.frame with only containing chromosomes in chrs.
#'
#' @keywords internal
#' @noRd
.chr_filter <- function(x, chrs) {

    # For R CMD check
    chr <- NULL

    # check if different chromosome formats
    if ((!any(chrs %in% x[["chr"]]))) {
        stop("No chromosomes in chrs match those of the junction data.")
    } else if (!(all(chrs %in% x[["chr"]]))) {
        warning(stringr::str_c(
            "The following chromosomes in chrs are not found in the junction data: ",
            stringr::str_c(chrs[!(chrs %in% x[["chr"]])], collapse = ", ")
        ))
    }

    x <- x %>% dplyr::filter(chr %in% chrs)

    return(x)
}

#' Load coverage for set of genomic regions
#'
#' `.coverage_load` loads coverage across a set of genomic regions from a
#' BigWig file using \href{https://github.com/ChristopherWilks/megadepth}{megadepth}
#'
#' @param coverage_path path to BigWig (or BAM - STILL UNTESTED) file.
#' @param regions [GRangesList-class][GenomicRanges::GRangesList-class] object.
#' @param chr_format NULL or one of "chr" or "no_chr", indicating the chromsome
#'   format used in the `coverage_path`. Will convert the `regions` to
#'   this format if they are not already.
#' @param sum_fun "mean", "sum", "max", "min" indicating the summary function to
#'   perform across the coverage for each region.
#' @param out_dir directory to save the the outputs for running `megadepth`.
#'
#' @return numeric vector of equal length to `regions` with each value
#'   corresponding to the coverage across one genomic interval in
#'   `regions`.
#'
#' @keywords internal
#' @noRd
.coverage_load <- function(coverage_path,
    regions,
    chr_format = NULL,
    sum_fun,
    out_file = tempfile(),
    method = "md") {

    ##### check user input #####

    if (!sum_fun %in% c("mean", "sum", "max", "min")) {
        stop("sum_fun must be one of: 'mean', 'sum', 'max' or 'min'")
    }

    ##### check chr format #####

    if (!is.null(chr_format)) {
        regions <- .chr_check(regions, chr_format)
    }

    ##### load coverage #####

    if (method == "md") {
        temp_regions_path <- stringr::str_c(out_file, ".bed")
        temp_coverage_prefix <- stringr::str_c(out_file, "_coverage")

        regions %>%
            as.data.frame() %>%
            dplyr::select(seqnames, start, end, strand) %>%
            dplyr::mutate(
                start = start - 1, # megadepth uses 0-based indexing
                dummy_1 = ".",
                dummy_2 = "."
            ) %>%
            readr::write_delim(temp_regions_path,
                delim = "\t",
                col_names = FALSE
            )

        coverage <- megadepth::get_coverage(
            bigwig_file = coverage_path,
            op = sum_fun,
            annotation = temp_regions_path,
            prefix = temp_coverage_prefix
        )

        coverage <- mcols(coverage)[["score"]]

        if (sum(coverage, na.rm = TRUE) == 0) {
            warning("Total AUC across all regions was 0. Make sure chromosome format matches between input regions and bigWig/BAM file.")
        }
    } else if (method == "rt") {
        coverage <- rtracklayer::import(coverage_path,
            which = regions,
            as = "NumericList"
        )
    }

    return(coverage)
}

#' Cache a file if it is not found locally
#'
#' `.file_cache` will use [BiocFileCache][BiocFileCache::BiocFileCache-class]
#' will cache the file for faster repeated retrival, if it not found locally
#' (i.e. a URL).
#'
#' @param file_path a path to file of interest.
#'
#' @return file_path of cached file or unchanged file_path if found locally.
#'
#' @keywords internal
#' @noRd
.file_cache <- function(file_path) {
    if (!file.exists(file_path)) {

        # suppress warning for tidyverse deprecated funs (select_() used over select()) in BiocFileCache
        # exact = TRUE means exact match required, if F then uses regex search
        suppressWarnings(
            file_path <- BiocFileCache::bfcrpath(BiocFileCache::BiocFileCache(ask = FALSE),
                file_path,
                exact = TRUE
            )
        )
    }

    return(file_path)
}

#' Splits a GRanges object by it's start and end.
#'
#' `.get_gr_for_start_end` takes a [GRanges][GenomicRanges::GRanges-class]
#' object and generates 2, one containing only the start co-ordinate and the
#' other, the end.
#'
#' @param gr a [GRanges][GenomicRanges::GRanges-class] object.
#'
#' @return list containing 2 [GRanges][GenomicRanges::GRanges-class] objects,
#'   each with 1 range corresponding to every range in the input. One contains
#'   start positions, the other ends.
#'
#' @keywords internal
#' @noRd
.get_start_end <- function(x) {
    x_start <- x
    end(x_start) <- start(x_start)

    x_end <- x
    start(x_end) <- end(x_end)

    x_start_end <- list(
        start = x_start,
        end = x_end
    )

    return(x_start_end)
}

# mark the midpoints to label with the count
.mark_mid <- function(x) {

    # the point with the lowest absolute diff from the mid
    mid_index <- which.min(abs(x - mean(x)))

    mid_marked <- logical(length = length(x))
    mid_marked[mid_index] <- TRUE

    return(mid_marked)
}

#' Merges [CharacterList][IRanges::AtomicList] objects
#'
#' `.merge_lists` merges two [CharacterList][IRanges::AtomicList] objects into
#' one through the element-wise concatenation of the vectors inside each list
#'
#' @param x a [CharacterList][IRanges::AtomicList] object.
#' @param y a [CharacterList][IRanges::AtomicList] object.
#'
#' @return [CharacterList][IRanges::AtomicList] object the identical length to x
#'   and y. The elements of this output being the concantenation of the
#'   correponding elements in x and y.
#'
#' @keywords internal
#' @noRd
.merge_CharacterList <- function(x, y) {
    if ((length(x) != length(y))) {
        stop("lengths of x and y should be identical!")
    }

    # set all groups to equal to the indexes of x/y
    all_groups <- seq_along(x) %>%
        as.character()
    names(x) <- all_groups
    names(y) <- all_groups

    # unlist x_y into a vector
    # obtain the names which are the original groups of each elements
    x_y <- c(x, y) %>% unlist()
    groups <- names(x_y)

    x_y_merged <-
        x_y %>%
        unname() %>%
        split(f = groups %>% factor(all_groups)) %>%
        CharacterList()

    return(x_y_merged)
}

#' Generate outlier scores
#'
#' `.outlier_score` will employ an isolation forest using `sklearn` through
#' [reticulate][reticulate::reticulate] to score junctions by how outlier-y they
#' look based on disruptions to junction count and associated regions of
#' coverage (or any other junction-level feature).
#'
#' @inheritParams junction_outlier_score
#' @param features data.frame with columns detailing features to be inputted
#'   into an outlier detection model.
#'
#' @return numeric vector containing outlier scores for each junction.
#'
#' @keywords internal
#' @noRd
.outlier_score <- function(features, ...) {
    cl <- basilisk::basiliskStart(env_sklearn)

    outlier_scores <- basilisk::basiliskRun(cl, function() {
        sklearn <- reticulate::import("sklearn")
        od_model <- sklearn$ensemble$IsolationForest()
        od_model <- od_model$set_params(...)
        od_model_params <- od_model$get_params()
        od_model <- od_model$fit(features)
        outlier_scores <- od_model$decision_function(features)

        suppressWarnings(
            print(stringr::str_c(
                Sys.time(), " - fitting outlier detection model with parameters: ",
                stringr::str_c(names(od_model_params), "=", unname(od_model_params)) %>%
                    stringr::str_c(collapse = ", ")
            ))
        )

        return(outlier_scores)
    })

    basilisk::basiliskStop(cl)

    return(outlier_scores)
}

#' Load reference annotation
#'
#' `.ref_load`` loads reference annotation using
#' [makeTxDbFromGFF][GenomicFeatures::makeTxDbFromGFF] if a character or leaves
#' `ref` unchanged if already a [TxDb-class][GenomicFeatures::TxDb-class].
#'
#' @inheritParams junction_annot
#'
#' @return a [TxDb-class][GenomicFeatures::TxDb-class] object.
#'
#' @keywords internal
#' @noRd
.ref_load <- function(ref) {
    if (!is(ref, "character") & !is(ref, "TxDb")) {
        stop("ref must either be a path to the .gtf/gff3 file or a pre-loaded TxDb object")
    }

    if (is(ref, "character")) {
        print(stringr::str_c(Sys.time(), " - Importing gtf/gff3 as a TxDb..."))

        # cache the gtf/gff3 for faster repeated retrieval
        ref <- .file_cache(ref)

        # import gtf using refGenome, needed to obtain the annotated splice junctions easily
        ref <- GenomicFeatures::makeTxDbFromGFF(ref)
    }

    return(ref)
}

#' Re-groups values in vector or list
#'
#' `.regroup` takes as input a vector or list and regroups the contents. Each
#' element in outputted list will correspond to a group. Values in `x` that fall
#' into the same group will be merged/concatenated together to form a vector in
#' the resulting list. Empty groups will be not be dropped and filled with
#' vectors of length = 0. Based on the function [split][S4Vectors::splitAsList].
#'
#' @param x an atomic vector or list (only tested on character and
#'   [CharacterList][IRanges::AtomicList]).
#' @param groups vector of the same length as `x` indicating which group each
#'   corresponding element of `x` is assigned.
#' @param all_groups vector with values all groups and the order that they
#'   should be.
#'
#' @return list of the same length and order as `all_groups`.
#'
#' @keywords internal
#' @noRd
.regroup <- function(x, groups, all_groups) {

    # set name of x to groups
    # so when unlisted will retain these names
    # if groups is an integer, will be converted to character through names()
    names(x) <- groups

    # unlist() to ensure list inputs are a vector
    # this won't change vector inputs
    x <- x %>% unlist()

    # reset the groups to ensure this matches the
    # length of the unlisted lists
    groups <- names(x)

    # split() will regroup the vector into a list
    # with each element of the list corresponding to a group
    # this list will be all_groups in length
    # with empty groups containing of a vector of length = 0
    x_regrouped <- x %>%
        unname() %>% # remove unecessary names from outputted list
        S4Vectors::split(
            f = groups %>% factor(all_groups),
            drop = FALSE
        )

    return(x_regrouped)
}

#' Calculate z-score for x from the distribution y
#'
#' `.zscore` calculates a z-score for each junction count from a patient sample
#' (`x`), indicating it's deviation from the distribution of a controls counts
#' (`y`) of the same junction.
#'
#' @param x numeric vector containing case counts/coverage for 1 junction.
#' @param y numeric vector containing control counts for 1 junction.
#' @param sd_const a numeric scalar to be added to all control standard
#'   deviations. This is to prevent infinate/NaN values occuring when the
#'   standard deviation of the control counts is 0.
#'
#' @return numeric vector of length equal to `x` containing the z-score for the
#'   case junctions.
#'
#' @keywords internal
#' @noRd
.zscore <- function(x, y, sd_const = 0.01) {
    x_score <- (x - mean(y)) / (sd(y) + sd_const)
    return(x_score)
}
