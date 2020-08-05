#' Annotate junctions using reference annotation
#'
#' \code{junction_annot} annotates junctions by 1. Whether their start and/or
#' end position precisely overlaps with an annotated exon boundary and 2.
#' Whether that junction matches any intron definition from existing annotation.
#' Using this information along with the strand, junctions are categorised into
#' "annotated", "novel_acceptor", "novel_donor", "novel_combo",
#' "novel_exon_skip", "ambig_gene" and "unannotated".
#'
#' @param junctions junction data as a
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   object.
#' @param ref either path to gtf/gff3 or object of class `TxDb` imported using
#'   \code{\link[GenomicFeatures]{makeTxDbFromGFF}}.
#'
#' @return junctions as a
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   object with additional `rowData` detailing overlapping
#'   genes/transcripts/exons and junction categories.
#'
#' @examples
#'
#' ref <-
#'     "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
#' if (!exists("junctions_annoted")) {
#'     junctions_annoted <-
#'         junction_annot(
#'             junctions_example,
#'             ref
#'         )
#' }
#' junctions_annoted
#' @export
junction_annot <- function(junctions, ref) {

    ##### Check user input is correct #####

    if (!(methods::isClass(junctions, "RangedSummarisedExperiment"))) {
        stop("junctions must be in a RangedSummarisedExperiment format")
    }

    ##### Extract annotated exons/junctions co-ordinates from gtf #####

    print(stringr::str_c(Sys.time(), " - Obtaining co-ordinates of annotated exons and junctions from gtf/gff3..."))

    ref <- .ref_load(ref)
    ref_exons <- ref %>% GenomicFeatures::exons(columns = c("gene_id", "tx_name", "exon_name"))
    ref_introns <- ref %>%
        GenomicFeatures::intronsByTranscript() %>%
        unlist()

    ##### Obtain annotation through overlapping introns/exons #####

    print(stringr::str_c(Sys.time(), " - Getting junction annotation using overlapping exons..."))

    junctions <- .junction_annot_ref(junctions, ref_introns, ref_exons)

    ##### Tidy annotation - collapse gene annotation to per junction and infer strand #####

    print(stringr::str_c(Sys.time(), " - Tidying junction annotation..."))

    junctions <- .junction_annot_tidy(junctions)

    ##### Derive junction categories using strand & overlapping exon annotation #####

    print(stringr::str_c(Sys.time(), " - Deriving junction categories..."))

    junctions <- .junction_cat(junctions)

    print(stringr::str_c(Sys.time(), " - done!"))

    return(junctions)
}

#' Extracts annotation from the reference gtf
#'
#' \code{.junction_ref_annot} will find whether each junctions start/end
#' precisely matches the end/start of annotated exons. Then, for each hit will
#' annotate the start/end of the junction with the strand/exon/transcript/gene
#' from the reference annotation.
#'
#' @inheritParams junction_annot
#'
#' @param ref_exons annotated exons obtained using
#'   \code\link{[GenomicFeatures](exons)}.
#' @param ref_introns annotated introns obtained using
#'   \code\link{[GenomicFeatures](intronsByTranscript)}.
#' @param ignore.strand used by \code\link{[GenomicRanges](findOverlaps)}.
#'
#' @return junctions with annotation.
#'
#' @keywords internal
#' @noRd
.junction_annot_ref <- function(junctions, ref_introns, ref_exons, ignore.strand = FALSE) {

    ##### Do junctions match introns from the reference annotation? #####

    junctions_intron_hits <- junctions %>%
        findOverlaps(ref_introns,
            type = "equal",
            ignore.strand = ignore.strand
        )

    mcols(junctions)[["in_ref"]] <- seq_along(junctions) %in% unique(queryHits(junctions_intron_hits))

    ##### Do junctions start/end overlap with an exon end/start? #####

    # match junctions to exon definitions
    start(junctions) <- start(junctions) - 1
    end(junctions) <- end(junctions) + 1

    # make a gr where each junction/exon is marked by only a start or end co-ordinate
    junctions_start_end <- .get_start_end(junctions)
    ref_exons_start_end <- .get_start_end(ref_exons)

    ref_col_names <- c(ref_exons %>% mcols() %>% colnames(), "strand", "exon_width")

    for (start_end in c("start", "end")) {

        # only get hits between junc start/exon end or junc end/exon start
        # the other way (e.g. junc end/exon end) should not happen (only 0.05% of the data)
        end_start <- ifelse(start_end == "start", "end", "start")

        junctions_exon_hits <- findOverlaps(
            query = junctions_start_end[[start_end]],
            subject = ref_exons_start_end[[end_start]],
            type = "equal",
            ignore.strand = ignore.strand
        )

        for (i in seq_along(ref_col_names)) {

            # extract the values to be used for annotation of junctions
            if (ref_col_names[i] == "strand") {
                ref_col_raw <- ref_exons %>%
                    strand() %>%
                    as.character()
            } else if (ref_col_names[i] == "exon_width") {
                ref_col_raw <- ref_exons %>%
                    width() %>%
                    as.integer()
            } else {
                ref_col_raw <- mcols(ref_exons)[[ref_col_names[i]]]
            }

            ref_col_tidy <-
                .regroup(
                    x = ref_col_raw[subjectHits(junctions_exon_hits)], # subset annotation by the
                    groups = queryHits(junctions_exon_hits), # group hits by the junction they overlap
                    all_groups = seq_along(junctions) # each junction is a group
                )

            if (is.character(ref_col_tidy[[1]])) {
                ref_col_tidy <- ref_col_tidy %>%
                    CharacterList() %>%
                    unique() # unique values for when junction start/end overlaps e.g. two exons
            } else {
                ref_col_tidy <- ref_col_tidy %>%
                    IRanges::IntegerList()
            }

            mcols(junctions)[[stringr::str_c(ref_col_names[i], "_", start_end)]] <- ref_col_tidy
        }
    }

    # convert junc co-ords back to intron definitions
    start(junctions) <- start(junctions) + 1
    end(junctions) <- end(junctions) - 1

    return(junctions)
}


#' Tidying junction annotation
#'
#' \code{.junction_annot_tidy} merges the gene and strand details from the start and
#' end into one column per junction. Then, combines strand information from the
#' original RNA-seq based and that from overlapping annotation.
#'
#' @inheritParams junction_annot
#'
#' @return junctions with tidy annotation.
#'
#' @keywords internal
#' @noRd
.junction_annot_tidy <- function(junctions, cols_to_merge = c("gene_id", "strand")) {

    ##### Collapse gene_id/strand annotation from start/end #####

    # collapse gene/strand columns to per junction
    # instead of per start/end for easier querying
    for (col in cols_to_merge) {
        mcols(junctions)[[stringr::str_c(col, "_junction")]] <-
            .merge_CharacterList(
                x = mcols(junctions)[[stringr::str_c(col, "_start")]],
                y = mcols(junctions)[[stringr::str_c(col, "_end")]]
            ) %>%
            unique()
    }

    ##### Tidy strand #####

    # replacing empty strands ("none") and those with >1 strand ("ambig_gene") with "*"
    # this ensures each vector in CharacterList is of length 1
    # so can be unlisted and length(chr_list) == length(unlist(chr_list))
    mcols(junctions)[["strand_junction"]][lengths(mcols(junctions)[["strand_junction"]]) == 0] <- "*"
    mcols(junctions)[["strand_junction"]][lengths(mcols(junctions)[["strand_junction"]]) > 1] <- "*"

    # compare the strand obtained from annotation strand to original strand
    # salvage situations when either original or annotation strand is "*" and the other is "+" or "-"
    strand_orig <- as.character(strand(junctions))
    strand_annot <- unlist(mcols(junctions)[["strand_junction"]])

    strand(junctions) <- dplyr::case_when(
        strand_orig == strand_annot ~ strand_orig,
        strand_orig == "*" & strand_annot != "*" ~ strand_annot,
        strand_annot == "*" & strand_orig != "*" ~ strand_orig,
        TRUE ~ NA_character_
    )

    if (any(is.na(strand(junctions)))) {
        stop("There should be no strands left as NA after tidying...")
    }

    # remove to avoid confusion between strand() and strand_junc
    mcols(junctions)[["strand_junc"]] <- NULL

    return(junctions)
}

#' Categorises junctions depending on reference annotation and strand
#'
#' \code{.junction_cat} categories junctions into "annotated", "novel_acceptor",
#' "novel_donor", "novel_combo", "novel_exon_skip", "ambig_gene" and "none"
#' using information from annotation and strand.
#'
#' @inheritParams junction_annot
#'
#' @return junctions with additional metadata detailling junction categories.
#'
#' @keywords internal
#' @noRd
.junction_cat <- function(junctions, ref_junc) {

    # store strand out for readability
    strand_junc <- as.character(strand(junctions))

    mcols(junctions)[["type"]] <-
        dplyr::case_when(
            mcols(junctions)[["in_ref"]] == TRUE ~ "annotated",
            lengths(mcols(junctions)[["gene_id_junction"]]) == 0 ~ "unannotated",
            lengths(mcols(junctions)[["gene_id_junction"]]) > 1 ~ "ambig_gene", # after these checks lengths(gene_id_junction) must equal 1
            lengths(mcols(junctions)[["gene_id_start"]]) > 0 & lengths(mcols(junctions)[["gene_id_end"]]) > 0 ~ "novel_combo",
            strand_junc == "+" & lengths(mcols(junctions)[["gene_id_start"]]) > 0 ~ "novel_acceptor",
            strand_junc == "-" & lengths(mcols(junctions)[["gene_id_start"]]) > 0 ~ "novel_donor",
            strand_junc == "+" & lengths(mcols(junctions)[["gene_id_end"]]) > 0 ~ "novel_donor",
            strand_junc == "-" & lengths(mcols(junctions)[["gene_id_end"]]) > 0 ~ "novel_acceptor",
            TRUE ~ NA_character_
        )

    if (any(is.na(mcols(junctions)[["type"]]))) {
        stop("There should be no junction categories left as NA after tidying...")
    }

    # split the novel_combo into novel_combo and novel_exon_skip
    # do this separately to save time - case_when() evaluates each condition across all junctions
    # since each column is called as mcols(junctions)[["col"]]
    # the below only checks for novel_exon_skip in the novel_combo subset
    # the two distinguished by whether the start/end of junction overlap a matching transcript
    mcols(junctions)[["index_tmp"]] <- seq_along(junctions)
    novel_combo <- junctions[mcols(junctions)[["type"]] == "novel_combo"]
    novel_exon_skip_indexes <- mcols(novel_combo)[["index_tmp"]][any(mcols(novel_combo)[["tx_name_start"]] %in% mcols(novel_combo)[["tx_name_end"]])]
    mcols(junctions)[["type"]][novel_exon_skip_indexes] <- "novel_exon_skip"
    mcols(junctions)[["index_tmp"]] <- NULL

    # set junction categories as factor with all possible levels
    mcols(junctions)[["type"]] <- mcols(junctions)[["type"]] %>%
        factor(levels = c(
            "annotated",
            "novel_acceptor",
            "novel_donor",
            "novel_exon_skip",
            "novel_combo",
            "ambig_gene",
            "unannotated"
        ))

    return(junctions)
}
