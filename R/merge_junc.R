#' Merge junction data from multiple samples
#'
#' \code{merge_junc} merges the junction data from several samples. Uses a
#' \code{\link[dplyr]{full_join}}, so retains all unique junctions from every
#' sample.
#'
#' @param junc_paths paths to files containing junction data.
#' @param sample_ids sample ids in the same order as \code{junc_paths}.
#' @param load_func function to load in junctions. By default, requires STAR
#'   formatted junctions. But this could be switched dependent on user-specific
#'   junction data format.
#' @param chr_to_filter chromosomes you would like to keep. By default, no
#'   filter is applied.
#'
#' @return junction co-ordinates and the raw counts from each sample.
#'
#' @export
merge_junc <- function(junc_paths, sample_ids, load_func = .load_STAR, chr_to_filter = NULL){

  if(length(junc_paths) != length(sample_ids)){

    stop("Number of junc_paths does not equal to the number of sample_ids")

  }

  for(i in seq_along(junc_paths)){

    print(stringr::str_c(Sys.time(), " - loading junctions for ", sample_ids[i], "..."))

    junc_df <- load_func(junc_paths[i], sample_ids[i])

    # filter for chromsomes
    if(!is.null(chr_to_filter)){

      # catch different chromosome formats
      if(any(chr_to_filter %in% junc_df[["chr"]]) == F){

        stop("No chromosomes in chr_to_filter match those of junction data")

      }

      junc_df <- junc_df %>%
        dplyr::filter(chr %in% chr_to_filter)

    }

    if(i == 1){

      junc_df_all <- junc_df

    }else{

      junc_df_all <- junc_df_all %>%
        dplyr::full_join(junc_df, by = c("chr", "start", "end", "strand"))

    }

  }

  # does the number unique chr_start_ends match the total length of junc_df?
  # if not, then warn user that some strands of junctions that have identical coords may be different
  uniq_junc_n <- junc_df_all %>%
    dplyr::group_by(chr, start, end) %>%
    dplyr::summarise(n = dplyr::n())

  if(any(uniq_junc_n[["n"]] != 1)){

    warning("Some juncs with identical co-ords occupy more than 1 row.\nPotentially due to junctions with idenitical co-ords, but differing strands.")

  }

  # replace all missing count values with 0
  junc_df_all[is.na(junc_df_all)] <- 0

  # structure as metadata and counts
  juncs <-
    list(metadata = junc_df_all %>%
           dplyr::select(chr, start, end, strand) %>%
           GRanges(),
         raw_count = junc_df_all %>%
           dplyr::select(-chr, -start, -end, -strand))

  print(stringr::str_c(Sys.time(), " - done!"))

  return(juncs)

}

# function used to load in junctions
# can be replaced if user juncs in a diff format
.load_STAR <- function(junc_path, sample_id){

  junc_df <-
    readr::read_delim(junc_path,
                      delim = "\t",
                      col_names = c("chr", "start", "end", "strand",
                                    "intron_motif", "annotation", "uniq_map_read_count", "multi_map_read_count", "max_overhang"),
                      col_types = readr::cols(chr = "c", .default = "i")) %>%
    dplyr::mutate(strand = dplyr::case_when(strand == 0 ~ "*",
                                            strand == 1 ~ "+",
                                            strand == 2 ~ "-")) %>%
    dplyr::select(chr:strand, !!sample_id := uniq_map_read_count)

  return(junc_df)

}
