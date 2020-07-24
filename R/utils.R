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
        GenomeInfoDb::seqlevels(x) <- GenomeInfoDb::seqlevels(x) %>%
            stringr::str_c("chr", .)
    } else if (chr_format == "no_chr" & any(stringr::str_detect(seqnames(x), "chr"))) {
        GenomeInfoDb::seqlevels(x) <- GenomeInfoDb::seqlevels(x) %>%
            stringr::str_replace("chr", "")
    }

    return(x)
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

    if (class(ref) == "character") {
        print(stringr::str_c(Sys.time(), " - Importing gtf/gff3 as a TxDb..."))

        # import gtf using refGenome, needed to obtain the annotated splice junctions easily
        ref <- GenomicFeatures::makeTxDbFromGFF(ref)
    }

    return(ref)
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

#' Re-groups values in vector or list
#'
#' \code{.regroup} takes as input a vector or list and regroups the
#' contents. Each element in outputted list will correspond to a group. Values
#' in \code{x} that fall into the same group will be merged/concatenated
#' together to form a vector in the resulting list. Empty groups will be not be
#' dropped and filled with vectors of length = 0. Based on the function
#' \code{\link{[S4Vectors](split)}}).
#'
#' @param x an atomic vector or list (only tested on character and
#'   \code{\link{[IRanges](CharacterList)}}).
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

#' Load coverage for set of genomic regions
#'
#' \code{.cov_load} loads coverage across a set of genomic regions from a BigWig
#' file using \url{https://github.com/ChristopherWilks/megadepth}
#'
#' @param regions [GRangesList-class][GenomicRanges::GRangesList-class] object.
#' @param cov_path path to BigWig (or BAM - STILL UNTESTED) file.
#' @param chr_format NULL or one of "chr" or "no_chr", indicating the chromsome
#'   format used in the \code{cov_path}. Will convert the \code{regions} to this
#'   format if they are not already.
#' @param sum_fun "mean", "sum", "max", "min" indicating the summary function to
#'   perform across the coverage for each region.
#'
#' @return numeric vector of equal length to \code{regions} with each value
#'   corresponding to the coverage across one genomic interval in
#'   \code{regions}.
#'
#' @keywords internal
#' @noRd
.cov_load <- function(regions, cov_path, chr_format = NULL, sum_fun) {

    ##### check user input #####

    if (!sum_fun %in% c("mean", "sum", "max", "min")) {
        stop("sum_fun must be one of: 'mean', 'sum', 'max' or 'min'")
    }

    ##### check chr format #####

    if (!is.null(chr_format)) {
        regions <- .chr_check(regions, chr_format)
    }

    ##### load coverage #####

    temp_regions_path <- stringr::str_c(tempdir(), "/regions.bed")
    temp_cov_prefix <- stringr::str_c(tempdir(), "/cov")

    regions %>%
        as.data.frame() %>%
        dplyr::select(seqnames, start, end) %>%
        dplyr::mutate(
            strand = "*",
            dummy_1 = ".",
            dummy_2 = "."
        ) %>%
        readr::write_delim(temp_regions_path, delim = "\t", col_names = FALSE)

    # check <- rtracklayer::import(cov_path, which = regions, as = "NumericList")

    system(
        command = stringr::str_c(
            "/tools/megadepth/megadepth ", cov_path,
            " --op ", sum_fun,
            " --annotation ", temp_regions_path,
            " ", temp_cov_prefix
        ),
        ignore.stdout = FALSE
    )

    suppressMessages(
        cov <- readr::read_delim(stringr::str_c(temp_cov_prefix, ".all.tsv"),
            delim = "\t",
            col_names = c("chr", "start", "end", "cov"),
            progress = FALSE
        )
    )

    cov <- cov[["cov"]]

    if (sum(cov, na.rm = TRUE) == 0) {
        warning("Total AUC across all regions was 0. Make sure chromsome format matches between input regions and bigWig/BAM file.")
    }

    return(cov)
}
