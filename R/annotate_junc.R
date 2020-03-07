#' Annotate junctions using existing annotation
#'
#' \code{annotate_junc_ref} annotates junctions by whether their start and end
#' position precisely overlap with a known exon boundary. On this basis,
#' categorises junctions into "annotated", "novel_acceptor", "novel_donor",
#' "novel_combo", "exon_skip", "ambig_gene" and "none".
#'
#' @param junc_metadata junction metadata output from \code{normalise_junc}. The
#'   essential component is the junction co-ordinates in a GRanges format.
#' @param gtf either path to gtf or of class \code{ensemblGenome} loaded using
#'   \code{\link{refGenome}}
#' @param ignore.strand
#'
#' @return list. Raw and normalised counts including metdata detailing junction
#'   annotation.
#'
#' @export
annotate_junc_ref <- function(junc_metadata, gtf, ignore.strand = F){

  ##### Extract annotated exons/junctions co-ordinates from gtf #####

  print(stringr::str_c(Sys.time(), " - Obtaining co-ords of annotated exons and junctions from gtf..."))

  if(class(gtf) == "character"){

    print(stringr::str_c(Sys.time(), " - Importing gtf..."))

    # import gtf using refGenome, needed to obtain the annotated splice junctions easily
    ref <- refGenome::ensemblGenome()
    refGenome::basedir(ref) <- dirname(gtf)
    refGenome::read.gtf(ref, gtf %>% stringr::str_replace("/.*/", ""))

  }else if(class(gtf) == "ensemblGenome"){

    ref <- gtf

  }

  ref_exons <- ref@ev$gtf[ref@ev$gtf$feature == "exon",] %>% GenomicRanges::GRanges()
  ref_junc <- refGenome::getSpliceTable(ref)
  ref_junc <- ref_junc@ev$gtf

  ##### Find overlaps between junctions and and known exon boundaries #####

  print(stringr::str_c(Sys.time(), " - Finding overlaps between junctions and annotated exons..."))

  # add 1 to end and subtract 1 from start of junc co-ords to match exon defs
  GenomicRanges::start(junc_metadata) <- GenomicRanges::start(junc_metadata) - 1
  GenomicRanges::end(junc_metadata) <- GenomicRanges::end(junc_metadata) + 1

  # make a gr where each junc/exon is only marked by it's start or end co-ordinate
  ref_exons_sep_start_end <- .get_gr_for_start_end(ref_exons)
  junc_coords_sep_start_end <- .get_gr_for_start_end(junc_metadata)

  # only get hits between junc start/exon end or junc end/exon start
  # the other way (e.g. junc end/exon end) should not happen (only 0.05% of the data)
  junc_start_exon_end_hits <- GenomicRanges::findOverlaps(query = junc_coords_sep_start_end$start,
                                                          subject = ref_exons_sep_start_end$end,
                                                          type = "equal",
                                                          ignore.strand = ignore.strand)

  junc_end_exon_start_hits <- GenomicRanges::findOverlaps(query = junc_coords_sep_start_end$end,
                                                          subject = ref_exons_sep_start_end$start,
                                                          type = "equal",
                                                          ignore.strand = ignore.strand)

  # convert junc co-ords back to intron definitions
  GenomicRanges::start(junc_metadata) <- GenomicRanges::start(junc_metadata) + 1
  GenomicRanges::end(junc_metadata) <- GenomicRanges::end(junc_metadata) - 1

  ##### Derive junction annotation/classification through the overlapping hits #####

  print(stringr::str_c(Sys.time(), " - Deriving junction annotation..."))

  junc_metadata <-
    .convert_hits_to_annot_ref(junc_metadata,
                               junc_exon_hits = junc_start_exon_end_hits,
                               ref_exons,
                               ref_cols = c("strand", "gene_name", "gene_id", "transcript_id", "exon_id"),
                               col_suffix = "_start")

  junc_metadata <-
    .convert_hits_to_annot_ref(junc_metadata,
                               junc_exon_hits = junc_end_exon_start_hits,
                               ref_exons,
                               ref_cols = c("strand", "gene_name", "gene_id", "transcript_id", "exon_id"),
                               col_suffix = "_end")

  print(stringr::str_c(Sys.time(), " - Tidying annotation..."))

  # obtain strand and gene combined per junction
  junc_metadata$strand_junc <- .merge_lists(junc_metadata$strand_start, junc_metadata$strand_end)
  junc_metadata$gene_name_junc <- .merge_lists(junc_metadata$gene_name_start, junc_metadata$gene_name_end)
  junc_metadata$gene_id_junc <- .merge_lists(junc_metadata$gene_id_start, junc_metadata$gene_id_end)

  # infer the annotation based on 1. the presence of junction in ref database and 2. using the strand acceptor/donor
  ref_junc_gr <- ref_junc %>% dplyr::rename(start = lend, end = rstart) %>% GenomicRanges::GRanges() %>% unique()
  start(ref_junc_gr) <- start(ref_junc_gr) + 1
  end(ref_junc_gr) <- end(ref_junc_gr) - 1

  suppressWarnings(annot_hits <- GenomicRanges::findOverlaps(query = junc_metadata, subject = ref_junc_gr, type = "equal"))

  junc_metadata$junc_in_ref <- 1:length(junc_metadata) %in% S4Vectors::queryHits(annot_hits)
  junc_metadata$annot_ref <- ifelse(lengths(junc_metadata$gene_name_junc) == 0, "none",
                                   ifelse(lengths(junc_metadata$gene_name_junc) >= 2, "ambig_gene",
                                          ifelse(all(junc_metadata$strand_junc == "+"),
                                                 ifelse(lengths(junc_metadata$strand_start) == 1,
                                                        ifelse(lengths(junc_metadata$strand_end) == 1,
                                                               ifelse(junc_metadata$junc_in_ref == T, "annotated", "exon_skip"), "novel_acceptor"), "novel_donor"),
                                                 ifelse(lengths(junc_metadata$strand_start) == 1,
                                                        ifelse(lengths(junc_metadata$strand_end) == 1,
                                                               ifelse(junc_metadata$junc_in_ref == T, "annotated", "exon_skip"), "novel_donor"), "novel_acceptor"))))

  # we can salvage some of those ambig_gene, juncs matching exons from 2 genes - ~13,000/800,000 junctions (1.6%)
  # by using the gene_info from the annotation
  annot_ambig_indexes <- which(junc_metadata$junc_in_ref == T & junc_metadata$annot_ref == "ambig_gene")
  annot_ambig_hits_df <- annot_hits %>% as.data.frame() %>% filter(queryHits %in% annot_ambig_indexes)

  junc_metadata$gene_name_junc <- replace(x = junc_metadata$gene_name_junc,
                                          list = annot_ambig_indexes,
                                          values = .convert_hits_to_list(x = ref_junc_gr$gene_name, y = annot_ambig_hits_df))

  junc_metadata$gene_id_junc <- replace(x = junc_metadata$gene_id_junc,
                                        list = annot_ambig_indexes,
                                        values = .convert_hits_to_list(x = ref_junc_gr$gene_id, y = annot_ambig_hits_df))

  junc_metadata$strand_junc <- replace(x = junc_metadata$strand_junc,
                                        list = annot_ambig_indexes,
                                        values = .convert_hits_to_list(x = GenomicRanges::strand(ref_junc_gr), y = annot_ambig_hits_df))

  junc_metadata$annot_ref[annot_ambig_indexes] <- "annotated"

  print(stringr::str_c(Sys.time(), " - done!"))

  return(junc_metadata)

}

#' Converts overlapping hits into annotation
#'
#' \code{.convert_hits_to_annot_ref}
#'
#' @inheritparams annotate_junc_ref
#' @param junc_exon_hits hits.
#' @param ref_exons gr.
#' @param ref_cols chr.
#' @param col_suffix chr.
#'
#' @return gr. Junction
.convert_hits_to_annot_ref <- function(junc_metadata, junc_exon_hits, ref_exons, ref_cols, col_suffix){

  if(any("strand" %in% ref_cols)){

    ref_exons$strand <- GenomicRanges::strand(ref_exons)

  }

  .get_chr_list_ref_col <- function(junc_exon_hits, ref_exons, ref_col){

    query_hits_fct <- S4Vectors::queryHits(junc_exon_hits) %>% factor(levels = 1:length(junc_metadata))

    junc_exon_hits_gr_list <-
      ref_exons %>%
      mcols() %>% # extract metadata as DataFrame
      .[[ref_col]] %>% # extract the col you want
      .[S4Vectors::subjectHits(junc_exon_hits)] %>% # subset the exons by the hits
      S4Vectors::split(query_hits_fct, drop = F) %>% # split into groups based on overlapping juncs
      CharacterList() %>%
      unique()

  }

  # get the genic properties from ref in a CharacterList form
  for(i in seq_along(ref_cols)){

    ref_col <- ref_cols[i]
    mcols(junc_metadata)[[str_c(ref_col, col_suffix)]] <- .get_chr_list_ref_col(junc_exon_hits, ref_exons, ref_col)

  }

  return(junc_metadata)

}

#' Merges two Characterlists into 1 with elementwise concatenation of the
#' vectors inside each list
#'
#' \code{.merge_lists}
#'
#' @param x list. Elements contain chr vectors.
#' @param y list. Elements contain chr vectors.
#'
#' @return list. Elements contain chr vectors concatenated between x and y
.merge_lists <- function(x, y){

  if(!identical(names(x), names(y))) stop("names of x and y lists should be identical!")

  x_y <- c(x, y) %>% unlist()

  x_y_merged <-
    x_y %>%
    unname() %>%
    S4Vectors::split(f = names(x_y) %>% factor(levels = names(x))) %>% # making this a factor keeps the order
    IRanges::CharacterList() %>%
    unique()

  return(x_y_merged)

}

.convert_hits_to_list <- function(x, y){

  x_list <-
    S4Vectors::split(x = x[y$subjectHits] %>% as.character(),
                     f = y$queryHits %>% factor(levels = y$queryHits %>% unique()))

  return(x_list)

}

#' Tests annotate_junc_ref for whether N juncs are matching the correct exons
#' and also whether the annotation for acceptor donor fits expections
#'
#' \code{test_annotate_junc_ref}
#'
#' @param junc_metadata_w_ref
#' @param gtf
#' @param n_junc_to_test
#'
#' @return
test_annotate_junc_ref <- function(junc_metadata_w_ref, gtf, n_junc_to_test = 100){

  if(class(gtf) == "character"){

    print(stringr::str_c(Sys.time(), " - Importing gtf..."))

    # import gtf using refGenome, needed to obtain the annotated splice junctions easily
    ref <- refGenome::ensemblGenome()
    refGenome::basedir(ref) <- dirname(gtf)
    refGenome::read.gtf(ref, gtf %>% stringr::str_replace("/.*/", ""))

  }else if(class(gtf) == "ensemblGenome"){

    ref <- gtf

  }

  ref_exons <- ref@ev$gtf[ref@ev$gtf$feature == "exon", ] %>% GenomicRanges::GRanges()
  ref_exons <- ref_exons %>% keepSeqlevels(GenomeInfoDb::seqlevels(junc_metadata_w_ref), pruning.mode = "coarse")

  junc_indexes_to_test <- sample(1:length(junc_metadata_w_ref), n_junc_to_test)

  identical_gene_check <- logical(length = n_junc_to_test)

  for(i in seq_along(junc_indexes_to_test)){

    junc_index_to_test <- junc_indexes_to_test[i]

    junc_to_test <- junc_metadata_w_ref[junc_index_to_test]

    exon_hits <-
      which(as.vector(GenomicRanges::seqnames(ref_exons) == GenomicRanges::seqnames(junc_to_test) &
                        GenomicRanges::end(ref_exons) == (GenomicRanges::start(junc_to_test) - 1)))

    ref_exons_hits <- ref_exons[exon_hits]

    identical_gene_check[i] <- identical((junc_to_test$gene_id_start %>% unlist() %>% unname() %>% sort()), (ref_exons_hits$gene_id %>% unique() %>% sort()))

  }

  annot_ref_none <- junc_metadata_w_ref[junc_metadata_w_ref$annot_ref == "none"]
  annot_ref_ambig <- junc_metadata_w_ref[junc_metadata_w_ref$annot_ref == "ambig"]
  annot_ref_annotated <- junc_metadata_w_ref[junc_metadata_w_ref$annot_ref == "annotated"]
  annot_ref_donor <- junc_metadata_w_ref[junc_metadata_w_ref$annot_ref == "donor"]
  in_ref_annotated <- junc_metadata_w_ref[junc_metadata_w_ref$junc_in_ref == 1]

  annot_ref_none_check <- length(annot_ref_none$gene_name_junc %>% unlist()) == 0
  annot_ref_ambig_check <- all(lengths(annot_ref_ambig$gene_id_junc) >= 2)
  annot_ref_annotated_check <- all((lengths(annot_ref_annotated$strand_start) != 0) & (lengths(annot_ref_annotated$strand_end) != 0))
  annot_ref_donor_check <- length(annot_ref_donor[all(annot_ref_donor$strand_junc == "-") & (lengths(annot_ref_donor$gene_id_start) != 0)]) == 0

  in_ref_annotated_check <- all((in_ref_annotated$annot_ref %>% unique()) == "annotated")

  return(c(all(identical_gene_check), annot_ref_none_check, annot_ref_ambig_check, annot_ref_annotated_check, annot_ref_donor_check))

}
