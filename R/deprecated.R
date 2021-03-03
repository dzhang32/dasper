#' Deprecated function: Annotate junctions using reference annotation
#'
#' `annotate_junc_ref` is deprecated, please use `junction_annot` instead.
#' annotates junctions by whether their start and end position precisely
#' overlaps with a known exon boundary. Using this information along with the
#' strand, junctions are categorised into "annotated", "novel_acceptor",
#' "novel_donor", "novel_combo", "novel_exon_skip", "ambig_gene" and "none".
#'
#' @param junc_metadata junction metadata in a
#'   [GRanges-class][GenomicRanges::GRanges-class] format, the essential
#'   component being the junction co-ordinates.
#' @param gtf either path to gtf or object of class \code{ensemblGenome} loaded
#'   using \code{refGenome}.
#'
#' @return junction metadata as a [GRanges-class][GenomicRanges::GRanges-class]
#'   object with additional columns that detail overlapping
#'   genes/transcripts/exons and junction categories.
#'
#' @keywords internal
#' @rdname deprecated
annotate_junc_ref <- function(junc_metadata, gtf) {
    lifecycle::deprecate_warn("0.99.0", "annotate_junc_ref()", "dasper::junction_annot()")

    ##### Check user input is correct #####

    if (!("GRanges" %in% class(junc_metadata))) stop("junction_metadata must be in a GRanges format")

    if (all(!(c("character", "ensemblGenome") %in% class(gtf)))) {
        stop("gtf must either be a path to the .gtf file or a pre-loaded gtf of class ensemblGenome")
    }

    ##### Extract annotated exons/junctions co-ordinates from gtf #####

    print(stringr::str_c(Sys.time(), " - Obtaining co-ordinates of annotated exons and junctions from gtf..."))

    if (class(gtf) == "character") {
        print(stringr::str_c(Sys.time(), " - Importing gtf..."))

        # import gtf using refGenome, needed to obtain the annotated splice junctions easily
        ref <- ensemblGenome()
        basedir(ref) <- dirname(gtf)
        read.gtf(ref, gtf %>% stringr::str_replace(".*/", ""))
    } else if (class(gtf) == "ensemblGenome") {
        ref <- gtf
    }

    ref_exons <- ref@ev$gtf[ref@ev$gtf$feature == "exon", ] %>% GRanges()
    ref_junc <- getSpliceTable(ref)
    ref_junc <- ref_junc@ev$gtf

    ##### Obtain annotation through overlapping exons #####

    print(stringr::str_c(Sys.time(), " - Getting junction annotation using overlapping exons..."))

    junc_metadata <- .get_ref_exons_annot(junc_metadata, ref_exons)

    ##### Tidy annotation - collapse gene annotation to per junction and infer strand #####

    print(stringr::str_c(Sys.time(), " - Tidying junction annotation..."))

    junc_metadata <- .tidy_junc_annot(junc_metadata)

    ##### Derive junction categories using strand & overlapping exon annotation #####

    print(stringr::str_c(Sys.time(), " - Deriving junction categories..."))

    junc_metadata <- .classify_junc(junc_metadata, ref_junc)

    print(stringr::str_c(Sys.time(), " - done!"))

    return(junc_metadata)
}

#' @keywords internal
#' @rdname deprecated
.get_ref_exons_annot <- function(
    junc_metadata,
    ref_exons,
    ref_cols = c("strand", "gene_name", "gene_id", "transcript_id", "exon_id")) {

    # match junctions to exon definitions
    start(junc_metadata) <- start(junc_metadata) - 1
    end(junc_metadata) <- end(junc_metadata) + 1

    # make a gr where each junc/exon is marked by only a start or end co-ordinate
    junc_start_end <- .get_gr_for_start_end(junc_metadata)
    ref_exons_start_end <- .get_gr_for_start_end(ref_exons)

    for (start_end in c("start", "end")) {

        # only get hits between junc start/exon end or junc end/exon start
        # the other way (e.g. junc end/exon end) should not happen (only 0.05% of the data)
        end_start <- ifelse(start_end == "start", "end", "start")

        # avoid seqlevel non-overlap warnings
        suppressWarnings(
            junc_exon_hits <- findOverlaps(
                query = junc_start_end[[start_end]],
                subject = ref_exons_start_end[[end_start]],
                type = "equal",
                ignore.strand = F
            )
        )

        # set junc_hits to factor with levels containing all junction indexes
        # so split(drop = F) keeps all junctions
        # not only those which precisely overlap an exon boundary
        junc_hits_fct <- queryHits(junc_exon_hits) %>%
            factor(levels = 1:length(junc_metadata))

        for (j in seq_along(ref_cols)) {

            # extract the values from exon metadata column of interest
            if (ref_cols[j] != "strand") {
                ref_col_values <- ref_exons %>%
                    mcols() %>%
                    .[[ref_cols[j]]]
            } else {

                # if strand extract strand
                ref_col_values <- strand(ref_exons)
            }

            mcols(junc_metadata)[stringr::str_c(ref_cols[j], "_", start_end)] <- ref_col_values %>%
                .[subjectHits(junc_exon_hits)] %>% # subset the exons by those that overlap juncs
                split(junc_hits_fct, drop = F) %>% # split into groups based on index of overlapping junc
                CharacterList() %>%
                unique() # parrallelised unique - remove duplicates when for example strand if junc overlaps >1 exon
        }
    }

    # convert junc co-ords back to intron definitions
    start(junc_metadata) <- start(junc_metadata) + 1
    end(junc_metadata) <- end(junc_metadata) - 1

    return(junc_metadata)
}

#' @keywords internal
#' @rdname deprecated
.tidy_junc_annot <- function(
    junc_metadata,
    cols_to_merge = c("strand", "gene_name", "gene_id")) {

    # collapse gene/strand columns to per junc instead of per start/end for easier querying
    for (col in cols_to_merge) {
        mcols(junc_metadata)[[stringr::str_c(col, "_junc")]] <-
            .merge_lists(
                mcols(junc_metadata)[[stringr::str_c(col, "_start")]],
                mcols(junc_metadata)[[stringr::str_c(col, "_end")]]
            )
    }

    # replacing empty strands ("none") and those with >1 strand ("ambig_gene") with "*"
    # this ensures each vector in CharacterList is of length 1
    # so can be unlisted and length(chr_list) == length(unlist(chr_list))
    mcols(junc_metadata)[["strand_junc"]][lengths(mcols(junc_metadata)[["strand_junc"]]) == 0] <- "*"
    mcols(junc_metadata)[["strand_junc"]][lengths(mcols(junc_metadata)[["strand_junc"]]) > 1] <- "*"
    strand_annot <- unlist(mcols(junc_metadata)[["strand_junc"]])

    # compare the strand obtained from annotation strand to original strand
    orig_strand <- as.character(strand(junc_metadata))

    # salvage situations when either original or annotation strand is "*" and the other is "+" or "-"
    strand(junc_metadata) <- dplyr::case_when(
        orig_strand == strand_annot ~ orig_strand,
        orig_strand == "*" & strand_annot != "*" ~ strand_annot,
        strand_annot == "*" & orig_strand != "*" ~ orig_strand
    )

    # remove to avoid confusion between strand() and strand_junc
    mcols(junc_metadata)[["strand_junc"]] <- NULL

    return(junc_metadata)
}

#' @keywords internal
#' @rdname deprecated
.classify_junc <- function(junc_metadata, ref_junc) {

    # find whether junction is found in splice table
    ref_junc_gr <- ref_junc %>%
        dplyr::rename(start = lend, end = rstart) %>%
        dplyr::mutate(
            start = start + 1, # match exon boundaries to intron co-ords
            end = end - 1
        ) %>%
        GRanges() %>%
        unique()

    # avoid diff seqlevels warning
    suppressWarnings(annot_hits <- findOverlaps(
        query = junc_metadata,
        subject = ref_junc_gr,
        type = "equal"
    ))

    mcols(junc_metadata)[["junc_in_ref"]] <- 1:length(junc_metadata) %in% queryHits(annot_hits)

    # classify junctions
    # separate strand out for readability
    strand_junc <- as.character(strand(junc_metadata))
    mcols(junc_metadata)[["junc_cat"]] <-
        dplyr::case_when(
            mcols(junc_metadata)[["junc_in_ref"]] == T ~ "annotated",
            lengths(mcols(junc_metadata)[["gene_name_junc"]]) == 0 ~ "none",
            lengths(mcols(junc_metadata)[["gene_name_junc"]]) > 1 ~ "ambig_gene", # after these checks lengths(gene_name_junc) must equal 1
            lengths(mcols(junc_metadata)[["gene_name_start"]]) > 0 & lengths(mcols(junc_metadata)[["gene_name_end"]]) > 0 ~ "novel_combo",
            strand_junc == "+" & lengths(mcols(junc_metadata)[["gene_name_start"]]) > 0 ~ "novel_acceptor",
            strand_junc == "-" & lengths(mcols(junc_metadata)[["gene_name_start"]]) > 0 ~ "novel_donor",
            strand_junc == "+" & lengths(mcols(junc_metadata)[["gene_name_end"]]) > 0 ~ "novel_donor",
            strand_junc == "-" & lengths(mcols(junc_metadata)[["gene_name_end"]]) > 0 ~ "novel_acceptor"
        )

    # split the novel_combo into novel_combo and novel_exon_skip
    # do this separately to save time - case_when evaluates each condition across all junctions
    # since each column is called as mcols(junc_metadata)[["col"]]
    # converting this data_frame() parses CharacterLists to list(), so also non-optimal solution
    mcols(junc_metadata)[["index_tmp"]] <- 1:length(junc_metadata)
    novel_combo <- junc_metadata[mcols(junc_metadata)[["junc_cat"]] == "novel_combo"]
    exon_skip_indexes <- mcols(novel_combo)[["index_tmp"]][any(novel_combo$transcript_id_start %in% novel_combo$transcript_id_end)]
    mcols(junc_metadata)[["junc_cat"]][exon_skip_indexes] <- "novel_exon_skip"
    mcols(junc_metadata)[["index_tmp"]] <- NULL

    return(junc_metadata)
}

#' @keywords internal
#' @rdname deprecated
.get_gr_for_start_end <- function(gr) {
    gr_start <- gr
    end(gr_start) <- start(gr_start)

    gr_end <- gr
    start(gr_end) <- end(gr_end)

    gr_start_end_list <- list(
        start = gr_start,
        end = gr_end
    )

    return(gr_start_end_list)
}

#' @keywords internal
#' @rdname deprecated
.merge_lists <- function(x, y) {
    if (!identical(names(x), names(y))) stop("names of x and y lists should be identical!")

    x_y <- c(x, y) %>% unlist()

    x_y_merged <-
        x_y %>%
        unname() %>%
        split(f = names(x_y) %>%
            factor(levels = names(x))) %>% # required to keep all levels/names
        CharacterList() %>%
        unique()

    return(x_y_merged)
}
