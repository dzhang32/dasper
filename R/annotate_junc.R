#' Annotate junctions using reference annotation
#'
#' \code{annotate_junc_ref} annotates junctions by whether their start and end
#' position precisely overlap with a known exon boundary. On this basis, then
#' categorises junctions into "annotated", "novel_acceptor", "novel_donor",
#' "novel_combo", "novel_exon_skip", "ambig_gene" and "none".
#'
#' @param junc_metadata junction metadata in a
#'   \code{\link[GenomicRanges]{GRanges}} format. The essential component is the
#'   junction co-ordinates. Other metadata columns will be returned unchanged.
#' @param gtf either path to gtf or object of class \code{ensemblGenome} loaded
#'   using \code{\link{refGenome}}.
#' @param ignore.strand whether to use the strand information when finding hits
#'   between exons and junctions.
#'
#' @return junction metadata with additional columns that detail overlapping
#'   genes/transcripts/exons and junction categories.
#'
#' @export
annotate_junc_ref <- function(junc_metadata, gtf, ignore.strand = F){

  ##### Check user input is correct #####

  if(!("GRanges" %in% class(junc_metadata))) stop("junction_metadata must be in a GRanges format")

  if(all(!(c("character", "ensemblGenome") %in% class(gtf)))){

    stop("gtf must either be a path to the .gtf file or a pre-loaded gtf of class ensemblGenome")

  }

  ##### Extract annotated exons/junctions co-ordinates from gtf #####

  print(stringr::str_c(Sys.time(), " - Obtaining co-ordinates of annotated exons and junctions from gtf..."))

  if(class(gtf) == "character"){

    print(stringr::str_c(Sys.time(), " - Importing gtf..."))

    # import gtf using refGenome, needed to obtain the annotated splice junctions easily
    ref <- refGenome::ensemblGenome()
    refGenome::basedir(ref) <- dirname(gtf)
    refGenome::read.gtf(ref, gtf %>% stringr::str_replace(".*/", ""))

  }else if(class(gtf) == "ensemblGenome"){

    ref <- gtf

  }

  ref_exons <- ref@ev$gtf[ref@ev$gtf$feature == "exon",] %>% GenomicRanges::GRanges()
  ref_junc <- refGenome::getSpliceTable(ref)
  ref_junc <- ref_junc@ev$gtf

  ##### Obtain junction annotation through overlapping exons #####

  print(stringr::str_c(Sys.time(), " - Getting junction annotation using overlapping exons..."))

  # obtain the reference annotation from the overlapping exons
  junc_metadata <- .get_ref_exons_annot(junc_metadata,
                                        ref_exons,
                                        ignore.strand = ignore.strand)

  ##### Tidy junction annotation #####

  print(stringr::str_c(Sys.time(), " - Tidying junction annotation..."))

  # collapse gene/strand columns to per junc instead of per start/end for easier querying
  for(col_to_merge in c("strand", "gene_name", "gene_id")){

    mcols(junc_metadata)[[stringr::str_c(col_to_merge, "_junc")]] <-
      .merge_lists(mcols(junc_metadata)[[stringr::str_c(col_to_merge, "_start")]],
                   mcols(junc_metadata)[[stringr::str_c(col_to_merge, "_end")]])

  }

  ##### Derive junction categories #####

  print(stringr::str_c(Sys.time(), " - Deriving junction categories..."))

  junc_metadata <- .classify_junc(junc_metadata, ref_junc)

  print(stringr::str_c(Sys.time(), " - done!"))

  return(junc_metadata)

}

#' Extracts annotation from the reference gtf
#'
#' \code{.convert_hits_to_annot_ref} will find the overlap of each junction with
#' the annotated exons. Then, for each corresponding hit, annotates each
#' junction with the strand/exon/transcript/gene from the reference annotation.
#'
#' @inheritParams annotate_junc_ref
#'
#' @param ref_exons annotated exons imported from the gtf.
#' @param junc_start_end start and end of the junctions returned from
#'   \code{\link{.get_gr_for_start_end}}.
#' @param ref_exons_start_end start and end of the exons returned from
#'   \code{\link{.get_gr_for_start_end}}.
#' @param ref_cols column names matching the annotation columns from the
#'   reference.
#'
#' @return junction metadata with annotation.
.get_ref_exons_annot <- function(junc_metadata, ref_exons, junc_start_end, ref_exons_start_end, ignore.strand,
                                 ref_cols = c("strand", "gene_name", "gene_id", "transcript_id", "exon_id")){

  # match junctions to exon definitions
  GenomicRanges::start(junc_metadata) <- GenomicRanges::start(junc_metadata) - 1
  GenomicRanges::end(junc_metadata) <- GenomicRanges::end(junc_metadata) + 1

  # make a gr where each junc/exon is marked by only a start or end co-ordinate
  junc_start_end <- .get_gr_for_start_end(junc_metadata)
  ref_exons_start_end <- .get_gr_for_start_end(ref_exons)

  for(start_end in c("start", "end")){

    # only get hits between junc start/exon end or junc end/exon start
    # the other way (e.g. junc end/exon end) should not happen (only 0.05% of the data)
    end_start <- ifelse(start_end == "start", "end", "start")

    # avoid seqlevel non-overlap warnings
    suppressWarnings(
      junc_exon_hits <- GenomicRanges::findOverlaps(query = junc_start_end[[start_end]],
                                                    subject = ref_exons_start_end[[end_start]],
                                                    type = "equal",
                                                    ignore.strand = ignore.strand)
    )

    # set junc_hits to factor with levels containing all junction indexes
    # so split(drop = F) keeps all junctions
    # not only those which precisely overlap an exon boundary
    junc_hits_fct <- queryHits(junc_exon_hits) %>%
      factor(levels = 1:length(junc_metadata))

    for(j in seq_along(ref_cols)){

      # extract the values from exon metadata column of interest
      if(ref_cols[j] != "strand"){

        ref_col_values <- ref_exons %>%
          mcols() %>%
          .[[ref_cols[j]]]

      }else{

        # if strand extract strand
        ref_col_values <- GenomicRanges::strand(ref_exons)

      }

      mcols(junc_metadata)[stringr::str_c(ref_cols[j], "_", start_end)] <- ref_col_values %>%
        .[subjectHits(junc_exon_hits)] %>% # subset the exons by those that overlap juncs
        split(junc_hits_fct, drop = F) %>% # split into groups based on index of overlapping junc
        IRanges::CharacterList() %>%
        unique() # parrallelised unique - remove duplicates in for example strand if junc overlaps >1 exon

    }
  }

  # convert junc co-ords back to intron definitions
  GenomicRanges::start(junc_metadata) <- GenomicRanges::start(junc_metadata) + 1
  GenomicRanges::end(junc_metadata) <- GenomicRanges::end(junc_metadata) - 1

  return(junc_metadata)

}

#' Classifies junctions
#'
#' \code{.classify_junc} categories junctions into "annotated",
#' "novel_acceptor", "novel_donor", "novel_combo", "exon_skip", "ambig_gene" and
#' "none" using information from annotation and strand. Adds two additional columns
#'
#' @inheritParams annotate_junc_ref
#'
#' @return
.classify_junc <- function(junc_metadata, ref_junc){

  # find whether junction is found in splice table
  ref_junc_gr <- ref_junc %>%
    dplyr::rename(start = lend, end = rstart) %>%
    dplyr::mutate(start = start + 1, # match exon boundaries to intron co-ords
                  end = end -1) %>%
    GenomicRanges::GRanges() %>%
    unique()

  # avoid diff seqlevels warning
  suppressWarnings(annot_hits <- GenomicRanges::findOverlaps(query = junc_metadata,
                                                             subject = ref_junc_gr,
                                                             type = "equal"))

  mcols(junc_metadata)["junc_in_ref"] <- 1:length(junc_metadata) %in% queryHits(annot_hits)

  # classify junctions
  mcols(junc_metadata)["junc_cat"] <-
    dplyr::case_when(junc_metadata$junc_in_ref == T ~ "annotated",
                     lengths(junc_metadata$gene_name_junc) == 0 ~ "none",
                     lengths(junc_metadata$gene_name_junc) > 1 ~ "ambig_gene", # after these checks lengths(gene_name_junc) must equal 1
                     lengths(junc_metadata$gene_name_start) > 0 & lengths(junc_metadata$gene_name_end) > 0 &
                       any(junc_metadata$transcript_id_start %in% junc_metadata$transcript_id_end) ~ "novel_exon_skip",
                     lengths(junc_metadata$gene_name_start) > 0 & lengths(junc_metadata$gene_name_end) > 0 ~ "novel_combo",
                     all(junc_metadata$strand_junc == "+") & lengths(junc_metadata$gene_name_start) > 0 ~ "novel_acceptor",
                     all(junc_metadata$strand_junc == "-") & lengths(junc_metadata$gene_name_start) > 0 ~ "novel_donor",
                     all(junc_metadata$strand_junc == "+") & lengths(junc_metadata$gene_name_end) > 0 ~ "novel_donor",
                     all(junc_metadata$strand_junc == "-") & lengths(junc_metadata$gene_name_end) > 0 ~ "novel_acceptor")

  return(junc_metadata)

}
