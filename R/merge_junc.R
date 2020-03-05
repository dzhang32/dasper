#' Merge juncs from multiple samples
#'
#' \code{merge_junc} merges the STAR SJ.out outputs from several samples into a
#' count df. Uses a full join, so keeps all unique junctions from every sample.
#'
#' @param junc_paths chr. Paths to the SJ.out files.
#' @param sample_ids chr. Sample ids in the same length and order and the
#'   junc_paths
#' @param chr_to_filter chr. Chromosomes you would like to keep.
#'
#' @return df. Junction co-ordinates and columns with the raw counts from each
#'   sample
#'
#' @export
merge_junc <- function(junc_paths, sample_ids, chr_to_filter = NULL){

  if(length(junc_paths) != length(sample_ids)){

    stop("The number of paths is not equal to the number of samples")

  }

  for(i in seq_along(junc_paths)){

    junc_path <- junc_paths[i]
    sample_id <- sample_ids[i]

    print(str_c(Sys.time(), " - loading junctions for ", sample_id, "..."))

    # remove strand info, STAR only uses intron motif to determine strand
    # instead, will add later when annotating by reference
    junc_df <-
      readr::read_delim(junc_path,
                        delim = "\t",
                        col_names = c("chr", "start", "end", "strand", "intron_motif", "annotation", "uniq_map_read_count", "multi_map_read_count", "max_overhang"),
                        col_types = cols(chr = "c", .default = "i")) %>%
      dplyr::select(chr:end, !!sample_id := uniq_map_read_count)

    if(!is.null(chr_to_filter)){

      if(!any(unique(junc_df$chr) %in% chr_to_filter)){

        stop("No junction chromosomes found in chr_to_filter")

      }

      junc_df <- junc_df %>% dplyr::filter(chr %in% chr_to_filter)

    }


    if(i == 1){

      junc_df_all <- junc_df

    }else{

      junc_df_all <-
        junc_df_all %>%
        dplyr::full_join(junc_df, by = c("chr", "start", "end"))

    }

  }

  print(str_c(Sys.time(), " - done!"))

  return(junc_df_all)

}
