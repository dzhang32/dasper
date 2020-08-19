#' Check chromosomes are correctly formatted
#'
#' \code{.chr_filter} will check whether the object \code{x} is matches the
#' format of chr_format. If not, will convert to the chromosome format.
#'
#' @param x object of [GRangesList-class][GenomicRanges::GRangesList-class] or
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   class.
#' @param chr_format one of "chr" or "no_chr".
#'
#' @return \code{x} with chromosomes matching \code{chr_format}.
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

#' Filter a df for chromosomes to keep
#'
#' \code{.chr_filter} will filter a df, keeping only chromosomes listed in the
#' chrs argument. It includes a couple of checks to notify users whether some or
#' all chrs are not found in their data. This aims to catch situations where
#' alternative nomenclature may be used - e.g. "M" vs "MT", or "chr1" vs "1".
#'
#' @param x df containing at minimum a column named "chr".
#' @inheritParams junction_load
#'
#' @return df with only containing chromosomes in chrs.
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
#' \code{.coverage_load} loads coverage across a set of genomic regions from a
#' BigWig file using \url{https://github.com/ChristopherWilks/megadepth}
#'
#' @param coverage_path path to BigWig (or BAM - STILL UNTESTED) file.
#' @param regions [GRangesList-class][GenomicRanges::GRangesList-class] object.
#' @param chr_format NULL or one of "chr" or "no_chr", indicating the chromsome
#'   format used in the \code{coverage_path}. Will convert the \code{regions} to
#'   this format if they are not already.
#' @param sum_fun "mean", "sum", "max", "min" indicating the summary function to
#'   perform across the coverage for each region.
#' @param out_dir directory to save the the outputs for running \code{megadepth}.
#'
#' @return numeric vector of equal length to \code{regions} with each value
#'   corresponding to the coverage across one genomic interval in
#'   \code{regions}.
#'
#' @keywords internal
#' @noRd
.coverage_load <- function(coverage_path, regions, chr_format = NULL, sum_fun, out_file = tempfile()) {

    ##### check user input #####

    if (!sum_fun %in% c("mean", "sum", "max", "min")) {
        stop("sum_fun must be one of: 'mean', 'sum', 'max' or 'min'")
    }

    ##### check chr format #####

    if (!is.null(chr_format)) {
        regions <- .chr_check(regions, chr_format)
    }

    ##### load coverage #####

    temp_regions_path <- stringr::str_c(out_file, ".bed")
    temp_coverage_prefix <- stringr::str_c(out_file, "_coverage")

    regions %>%
        as.data.frame() %>%
        dplyr::select(seqnames, start, end, strand) %>%
        dplyr::mutate(
            start = start - 1, # megadepth uses python indexing, -1 here to match rtracklayer
            dummy_1 = ".",
            dummy_2 = "."
        ) %>%
        readr::write_delim(temp_regions_path, delim = "\t", col_names = FALSE)

    # check <- rtracklayer::import(coverage_path, which = regions, as = "NumericList")

    system(
        command = stringr::str_c(
            "megadepth ", coverage_path,
            " --op ", sum_fun,
            " --annotation ", temp_regions_path,
            " ", temp_coverage_prefix
        ), ignore.stdout = TRUE
    )

    suppressMessages(
        coverage <- readr::read_delim(stringr::str_c(temp_coverage_prefix, ".all.tsv"),
            delim = "\t",
            col_names = c("chr", "start", "end", "cov"),
            col_types = readr::cols(chr = "c", start = "i", end = "i", cov = "n"),
            progress = FALSE
        )
    )

    coverage <- coverage[["cov"]]

    if (sum(coverage, na.rm = TRUE) == 0) {
        warning("Total AUC across all regions was 0. Make sure chromsome format matches between input regions and bigWig/BAM file.")
    }

    return(coverage)
}

#' Cache a file if it is not found locally
#'
#' \code{.file_cache} \code{\link[BiocFileCache](BiocFileCache)} will cache the
#' file for faster repeated retrival, if it not found locally (i.e. a URL).
#'
#' @param file_path a path to file of interest.
#'
#' @return \code{file_path} of cached file or unchanged \code{file_path} if
#'   found locally.
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

#' Splits a \code{GRanges} object by it's start and end.
#'
#' \code{.get_gr_for_start_end} takes a \code{GRanges} object and generates 2,
#' one containing only the start co-ordinate and the other, the end.
#'
#' @param gr Any \code{GRanges} object.
#'
#' @return list of 2 grs, each with 1 range corresponding to every range in the
#'   input. One contains start positions, the other ends.
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

#' Merges two \code{\link{[IRanges](CharacterList)}}s into one through the
#' element-wise concatenation of the vectors inside each list
#'
#' \code{.merge_lists}
#'
#' @param x \code{\link{[IRanges](CharacterList)}}
#' @param y \code{\link{[IRanges](CharacterList)}}
#'
#' @return \code{\link{[IRanges](CharacterList)}} the identical length to x and
#'   y. The elements of this output being the concantenation (\code{c()}) of the
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
#' \code{.outlier_score} will use an isolation forest to score junctions by how
#' abnormal they look based on disruptions to junction count and associated
#' regions of coverage.
#'
#' @inheritParams junction_outlier_score
#' @param features data.frame with columns detailing features to be inputted
#'   into an outlier detection model.
#'
#' @return outlier scores for each junction.
#'
#' @keywords internal
#' @noRd
.outlier_score <- function(features, ...) {

    # sklearn should be loaded in using .onLoad()
    od_model <- sklearn$ensemble$IsolationForest()

    # set parameters of
    od_model <- od_model$set_params(...)
    od_model_params <- od_model$get_params()

    suppressWarnings(
        print(stringr::str_c(
            Sys.time(), " - fitting outlier detection model with parameters: ",
            stringr::str_c(names(od_model_params), "=", unname(od_model_params)) %>%
                stringr::str_c(collapse = ", ")
        ))
    )

    od_model <- od_model$fit(features)
    outlier_scores <- od_model$decision_function(features)

    return(outlier_scores)
}

#' Load reference annotation
#'
#' \code{.ref_load} loads reference annotation using
#' \code{\link[GenomicFeatures]{makeTxDbFromGFF}} if a character or leaves
#' \code{ref} unchanged if already a [TxDb-class][GenomicFeatures::TxDb-class].
#'
#' @inheritParams junction_annot
#'
#' @return [TxDb-class][GenomicFeatures::TxDb-class] object.
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
#' \code{.regroup} takes as input a vector or list and regroups the
#' contents. Each element in outputted list will correspond to a group. Values
#' in \code{x} that fall into the same group will be merged/concatenated
#' together to form a vector in the resulting list. Empty groups will be not be
#' dropped and filled with vectors of length = 0. Based on the function
#' \code{\link[S4Vectors](split)}).
#'
#' @param x an atomic vector or list (only tested on character and
#'   \code{\link[IRanges](CharacterList)}).
#' @param groups vector of the same length as \code{x} indicating which group
#'   each corresponding element of \code{x} is assigned.
#' @param all_groups vector with values all groups and the order that they
#'   should be.
#'
#' @return list of the same length and order as \code{all_groups}.
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
#' \code{.zscore} calculates a z-score for each junction count from a patient
#' sample (\code{x}), indicating it's deviation from the distribution of a
#' controls counts (\code{y}) of the same junction.
#'
#' @param x numeric vector containing case counts/coverage for 1 junction.
#' @param y numeric vector containing control counts for 1 junction.
#' @param sd_const a numeric scalar to be added to all control standard
#'   deviations. This is to prevent infinate/NaN values occuring when the
#'   standard deviation of the control counts is 0.
#'
#' @return numeric vector of length equal to \code{x} containing the z-score for
#'   the case junctions.
#'
#' @keywords internal
#' @noRd
.zscore <- function(x, y, sd_const = 0.01) {
    x_score <- (x - mean(y)) / (sd(y) + sd_const)
    return(x_score)
}
