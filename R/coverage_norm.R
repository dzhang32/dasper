#' For each junction, obtain the normalised coverage across associated
#' exonic/intronic regions
#'
#' \code{coverage_norm} obtains regions of interest for each junction where
#' coverage disruptions would be expected. If ends of junctions are annotated,
#' it will use the overlapping exon definitions, picking the shortest exon when
#' multiple overlap one end. If unannotated, will use a user-defined width (i.e.
#' 20bp). To compare between samples (case/controls) coverage is normalised to a
#' set region. By default, the boundaries of each gene associated to a junction
#' are used as the region to normalise to.
#'
#' @inheritParams junction_annot
#'
#' @param unannot_width integer scalar determining the width of the region to
#'   obtain coverage from when the end of of a junction does not overlap an
#'   existing exon.
#' @param coverage_paths_case paths to the BigWig/BAM files containing the
#'   coverage of your case samples. Must be the same length and order to the
#'   rows of the assays in \code{junctions}.
#' @param coverage_paths_control paths to the BigWig/BAM files
#' @param coverage_chr_control either "chr" or "no_chr", indicating the
#'   chromosome format of control coverage data. To be used if you know the
#'   chromosome format of the control BigWig/BAM files is different to that of
#'   your junctions.
#'
#' @return list containing sublists, one for cases and the other controls. Each
#'   sublist contains 3 matrices, corresponding the coverage for each sample
#'   across the 2 exons and intron associated with every junction.
#'
#' @examples
#'
#' \dontrun{
#' # leave this as not run for now to save time for R CMD check
#' ref <- "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
#' ref <- GenomicFeatures::makeTxDbFromGFF(ref)
#' coverage_paths_case <- list.files("/data/RNA_seq_diag/mito/bw/", full.names = T)[1:2]
#' coverage_paths_control <-
#'     list.files(
#'         "/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/",
#'         full.names = T
#'     )[1:2]
#' coverage <- coverage_norm(
#'     junctions_annot_example,
#'     ref,
#'     unannot_width = 20,
#'     coverage_paths_case,
#'     coverage_paths_control,
#'     coverage_chr_control = "chr"
#' )
#' junctions
#' }
#'
#' @export
coverage_norm <- function(junctions, ref, unannot_width = 20, coverage_paths_case, coverage_paths_control, coverage_chr_control = NULL) {

    ##### Check user-input #####

    if (sum(colData(junctions)[["case_control"]] == "case") != length(coverage_paths_case)) {
        stop("Number of cases must equal the length of coverage_paths_case")
    }

    if (length(coverage_paths_control) < 2) {
        stop("coverage_paths_control must coverageer at least 2 controls")
    }

    if (!(coverage_chr_control %in% c("chr", "no_chr"))) {
        stop("coverage_chr_control must be one of 'chr' or 'no_chr'")
    }

    ##### Get exon/intronic regions for each junction #####

    print(stringr::str_c(Sys.time(), " - Obtaining exonic and intronic regions to load coverage from..."))

    coverage_regions <- .coverage_exon_intron(junctions, unannot_width)

    ##### Get region to normalise coverage with respect to #####

    print(stringr::str_c(Sys.time(), " - Obtaining regions to use to normalise coverage..."))

    # load in reference annotation
    ref <- .ref_load(ref)

    coverage_regions <- .coverage_norm_region(junctions, ref, coverage_regions)

    ##### Load coverage #####

    case_control_coverage <- .coverage_case_control_load(
        coverage_regions,
        coverage_paths_case,
        coverage_paths_control,
        coverage_chr_control
    )

    ##### Normalise coverage #####

    print(stringr::str_c(Sys.time(), " - Normalising coverage..."))

    case_control_coverage <- .coverage_norm(case_control_coverage)

    print(stringr::str_c(Sys.time(), " - done!"))

    return(case_control_coverage)
}

#' Obtain the co-ordinates corresponding to exons and introns for each junction
#'
#' \code{.coverage_exon_intron} obtains the exonic and intronic regions
#' corresponding each junction. If the junctions ends are annotated, it just
#' takes the exon definitions. If unannotated, then will use an arbitrary exon
#' width defined by the user in \code{unannot_width}. If multiple exons overlap
#' a junction end, the shortest definition is used. Intron definitions are the
#' junction co-ordinates.
#'
#' @inheritParams coverage_norm
#'
#' @return [GRangesList-class][GenomicRanges::GRangesList-class] containing 3
#'   sets of ranges, corresponding the the flanking exons and intron
#'   definitions.
#'
#' @keywords internal
#' @noRd
.coverage_exon_intron <- function(junctions, unannot_width) {
    coverage_regions <- GenomicRanges::GenomicRangesList()

    ##### Get exonic regions of interest #####

    for (start_end in c("start", "end")) {
        exon_widths <- mcols(junctions)[[stringr::str_c("exon_width_", start_end)]]

        # fill unannotated ends with user-defined width
        exon_widths <- replace(
            x = exon_widths,
            which(lengths(exon_widths) == 0),
            values = unannot_width
        ) %>%
            IRanges::IntegerList()

        stopifnot(all(lengths(exon_widths) > 0))

        # take the smallest exon if multiple overlap
        exon_widths <- min(exon_widths)

        # take off 1bp from the exon widths so resultant exon width remains the same
        if (start_end == "start") {
            exon_ends <- start(junctions) - 1
            exon_starts <- exon_ends - (exon_widths - 1)
        } else {
            exon_starts <- end(junctions) + 1
            exon_ends <- exon_starts + (exon_widths - 1)
        }

        coverage_regions[[stringr::str_c("exon_coords_", start_end)]] <-
            GRanges(
                seqnames = seqnames(junctions),
                ranges = IRanges::IRanges(
                    start = exon_starts,
                    end = exon_ends
                ),
                strand = strand(junctions)
            )
    }

    ##### Add intron/junction coords #####

    coverage_regions[["intron_coords"]] <-
        SummarizedExperiment::rowRanges(junctions)

    mcols(coverage_regions[["intron_coords"]]) <- NULL

    return(coverage_regions)
}

#' Obtain the region to normalise coverage
#'
#' \code{.coverage_norm_region} obtains the region to use to normalise coverage
#' across exons and intron for each junction. Currently, this is the gene
#' definition associated with each junction. For junctions that are unannotated,
#' the largest local region is taken from the start of the upstream to the end
#' of downstream exon. For ambiguous genes (> 1 gene associated to a junction),
#' the shortest gene definition is used.
#'
#' @inheritParams coverage_norm
#'
#' @param coverage_regions
#'
#' @return [GRangesList-class][GenomicRanges::GRangesList-class] containing an
#'   additional set of ranges containing regions to use to normalise.
#'
#' @keywords internal
#' @noRd
.coverage_norm_region <- function(junctions, ref, coverage_regions) {
    ref_genes <- ref %>%
        GenomicFeatures::genes(columns = c("gene_id", "exon_name")) %>%
        as.data.frame()

    coverage_norm_regions <- GenomicRanges::GenomicRangesList()

    ##### no gene #####

    no_gene_indexes <- which(lengths(mcols(junctions)[["gene_id_junction"]]) == 0)

    # for these we obtain the largest local region possible
    # from start exon overlapping start to end of exon overlapping end
    no_gene_coords <-
        GRanges(
            seqnames = seqnames(junctions)[no_gene_indexes],
            ranges = IRanges::IRanges(
                start = coverage_regions[["exon_coords_start"]][no_gene_indexes] %>% start(),
                end = coverage_regions[["exon_coords_end"]][no_gene_indexes] %>% end()
            ),
            strand = strand(junctions)[no_gene_indexes]
        )

    mcols(no_gene_coords)[["index"]] <- no_gene_indexes

    ##### single gene #####

    single_gene_indexes <- which(lengths(mcols(junctions)[["gene_id_junction"]]) == 1)

    # easiest case, we just use the gene coords
    single_gene_coords <- dplyr::tibble(
        gene_id = unlist(mcols(junctions)[["gene_id_junction"]][single_gene_indexes]),
        index = single_gene_indexes
    ) %>%
        dplyr::left_join(ref_genes, by = "gene_id") %>%
        dplyr::select(seqnames, start, end, strand, index) %>%
        GRanges()

    ##### ambig gene #####

    ambig_gene_indexes <- which(lengths(mcols(junctions)[["gene_id_junction"]]) > 1)

    ambig_gene_coords <- mcols(junctions)[["gene_id_junction"]][ambig_gene_indexes]
    names(ambig_gene_coords) <- ambig_gene_indexes

    # right now, this takes the smallest gene
    # however could be updated to be more complex
    ambig_gene_coords <- dplyr::tibble(
        index = names(unlist(ambig_gene_coords)),
        gene_id = unlist(ambig_gene_coords)
    ) %>%
        dplyr::left_join(ref_genes, by = "gene_id") %>%
        dplyr::group_by(index) %>%
        dplyr::filter(width == min(width)) %>%
        dplyr::select(seqnames, start, end, strand, index) %>%
        GRanges()

    ##### merge all gene coords together ####

    norm_coords <- c(
        no_gene_coords,
        single_gene_coords,
        ambig_gene_coords
    )

    # arrange the coords by index so they're order matches junctions
    mcols(norm_coords)[["index"]] <- mcols(norm_coords)[["index"]] %>% as.integer()
    norm_coords <- norm_coords %>% plyranges::arrange(index)
    stopifnot(identical(mcols(norm_coords)[["index"]], seq_along(junctions)))
    mcols(norm_coords)[["index"]] <- NULL

    coverage_regions[["norm_coords"]] <- norm_coords

    return(coverage_regions)
}

#' Loads coverage from case and control BigWig/BAM files
#'
#' \code{.coverage_case_control_load} uses \code{dasper:::.coverage_load} to
#' load coverage from BigWig/BAM files. It does this for both case and controls
#' and for each set of ranges in \code{coverage_regions}.
#'
#' @inheritParams coverage_norm
#' @inheritParams .coverage_norm_region
#'
#' @return list containing two lists called "case" and "control", each
#'   containing coverage matrices. The number of matrices is equal to the length
#'   of \code{coverage_regions}.
#'
#' @keywords internal
#' @noRd
.coverage_case_control_load <- function(coverage_regions, coverage_paths_case, coverage_paths_control, coverage_chr_control) {
    case_control_coverage <- list()

    for (case_control in c("case", "control")) {
        if (case_control == "case") {
            coverage_paths <- coverage_paths_case
            chr_format <- NULL # assume case coverage paths always same chr as junctions
        } else {
            coverage_paths <- coverage_paths_control
            chr_format <- coverage_chr_control
        }

        coverage_mats <- list()

        for (i in seq_along(coverage_regions)) {

            # sum to obtain total AUC for normalisation
            # mean for exon/intron regions
            sum_fun <- ifelse(names(coverage_regions)[i] == "norm_coords", "sum", "mean")

            coverage_mat <-
                matrix(
                    nrow = length(coverage_regions[[i]]),
                    ncol = length(coverage_paths)
                )

            for (j in seq_along(coverage_paths)) {
                print(stringr::str_c(
                    Sys.time(), " - Loading coverage across ", names(coverage_regions)[i],
                    " for ", case_control, " ", j, "/", length(coverage_paths), "..."
                ))

                coverage_mat[, j] <- .coverage_load(
                    regions = coverage_regions[[i]],
                    coverage_path = coverage_paths[j],
                    chr_format = chr_format,
                    sum_fun = sum_fun
                )
            }

            coverage_mats[[i]] <- coverage_mat
        }

        case_control_coverage[[case_control]] <- coverage_mats

        # convert names from coords to coverage to represent contents
        names(case_control_coverage[[case_control]]) <- names(coverage_regions) %>%
            stringr::str_replace("coords", "coverage")
    }

    return(case_control_coverage)
}

#' Normalises coverage across the exons/intron associated with each junction
#'
#' \code{.coverage_norm} normalises coverage across exons/introns by dividing
#' their coverage by the total AUC across the norm region.
#'
#' @param case_control_coverage list containing coverage for case and controls
#'   returned by \link{.coverage_case_control_load}
#' @param norm_const integer to add to the normalisation coverages. This
#'   prevents dividing by 0 and NaN/Inf values resulting.
#'
#' @return list containing two lists called "case" and "control". Each
#'   containing 3 matrices with normalised coverage across exons/introns for
#'   each junction/sample.
#'
#' @keywords internal
#' @noRd
.coverage_norm <- function(case_control_coverage, norm_const = 1) {

    # loop across case and controls
    for (case_control in c("case", "control")) {

        # and each coverage matrix within case/controls
        for (j in seq_along(case_control_coverage[[case_control]])) {

            # skip normalisation coverage as this does not need to be normalised itself
            if (names(case_control_coverage[[case_control]][j]) == "norm_coverage") {
                next
            } else {

                # normalise coverage across exons/intron
                # by dividing by the coverage across normalisation regions
                case_control_coverage[[case_control]][[j]] <-
                    case_control_coverage[[case_control]][[j]] / (case_control_coverage[[case_control]][["norm_coverage"]] + norm_const)

                case_control_coverage[[case_control]][[j]][is.na(case_control_coverage[[case_control]][[j]])] <- 0
            }
        }

        # remove norm coverage
        case_control_coverage[[case_control]] <-
            case_control_coverage[[case_control]][names(case_control_coverage[[case_control]]) != "norm_coverage"]
    }

    return(case_control_coverage)
}