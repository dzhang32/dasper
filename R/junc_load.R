#' Load junctions from patient and control RNA-seq data
#'
#' \code{junc_load} loads in raw patient and control junction data and formats
#' it into a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' Control samples can be user-inputted or selected from GTEx data publicly
#' released through the recount2 project
#' (\url{https://jhubiostatistics.shinyapps.io/recount/}) and downloaded
#' through snaptron (\url{http://snaptron.cs.jhu.edu/}). By default, \code{junc_load}
#' expects the junction data to be in STAR aligned format (SJ.out).
#'
#' @param junc_paths file path(s) to junction data.
#' @param metadata dataframe containing sample metadata with rows in the same
#'   order and corresponding to file path(s). Will be used as the \code{colData}
#'   of the \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param controls either a logical vector of the same length as paths with
#'   TRUEs labelling controls Or, "fibroblasts" representing the samples of
#'   which GTEx tissue to use as controls. If left unchanged, by default will
#'   assume all samples are patients.
#' @param load_func function to load in junctions. By default, requires STAR
#'   formatted junctions (SJ.out). But this can be switched dependent on the
#'   format of the user's junction data. Function must take as input a junction
#'   path then return a dataframe with the columns "chr", "start", "end",
#'   "strand" and "count".
#'
#' @inheritParams .chr_filter
#'
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing junction data.
#'
#' @examples
#'
#' \dontrun{
#'
#' example_juncs_1_path <-
#'     system.file("extdata", "example_juncs_1.txt",
#'         "dasper",
#'         mustWork = TRUE
#'     )
#'
#' juncs <-
#'     junc_load(
#'         junc_paths = c(example_juncs_1_path),
#'         metadata = dplyr::tibble(samp_id = c("example_1"))
#'     )
#'
#' juncs
#'
#' }
#'
#' @export
junc_load <- function(junc_paths,
    metadata = NULL,
    controls = rep(FALSE, length(junc_paths)),
    load_func = .load_STAR,
    chrs = NULL) {

    ##### Read in and merge junction data #####

    junc_df_all <- dplyr::tibble()

    for (i in seq_along(junc_paths)) {
        print(stringr::str_c(Sys.time(), " - loading junctions for sample ", i, "/", length(junc_paths), "..."))

        junc_df <- load_func(junc_paths[i])

        if (!is.null(chrs)) {
            junc_df <- .chr_filter(junc_df, chrs)
        }

        junc_df_all <- .junc_merge(junc_df_all, junc_df, i)
    }

    ##### Add control data/identifier #####

    if (any(controls %in% c("fibroblasts", "whole_blood"))) {



    } else {
        if (!(is.logical(controls)) | length(controls) != length(junc_paths)) {
            stop("Controls argument must be a logical vector of same length as the number of junction paths")
        }

        metadata <- metadata %>%
            dplyr::mutate(case_control = ifelse(controls, "control", "case"))
    }

    ##### Tidy junction data #####

    # replace all missing count (NA) values with 0
    junc_df_all[is.na(junc_df_all)] <- 0

    # convert junctions into a RangedSummarizedExperiment
    raw_counts <- junc_df_all %>%
        dplyr::select(-chr, -start, -end, -strand) %>%
        as.matrix()

    junc_coords <- junc_df_all %>%
        dplyr::select(chr, start, end, strand) %>%
        GRanges()

    juncs <-
        SummarizedExperiment::SummarizedExperiment(
            assays = list(raw_counts),
            rowRanges = junc_coords,
            colData = metadata
        )

    print(stringr::str_c(Sys.time(), " - done!"))

    return(juncs)
}

#' Load raw junction data
#'
#' \code{.load_STAR} will load raw junction data that is outputted from STAR
#' (SJ.out) into R. This will format the junction data, retaining only chr,
#' start, end, strand and count (uniq_map_read_count) columns.
#'
#' @param junc_path path to the junction data.
#'
#' @return df detailing junction co-ordinates and counts.
.load_STAR <- function(junc_path) {
    junc_df <-
        readr::read_delim(junc_path,
            delim = "\t",
            col_names = c(
                "chr", "start", "end", "strand",
                "intron_motif", "annotation", "uniq_map_read_count", "multi_map_read_count", "max_overhang"
            ),
            col_types = readr::cols(chr = "c", .default = "i")
        ) %>%
        dplyr::mutate(strand = dplyr::case_when(
            strand == 0 ~ "*",
            strand == 1 ~ "+",
            strand == 2 ~ "-"
        )) %>%
        dplyr::select(chr:strand, count := uniq_map_read_count)

    return(junc_df)
}

#' Merge two junction datasets together
#'
#' \code{.junc_merge} will merge two sets of junction data together. It uses a
#' full_join so will keep all rows from both datasets. It will also ensure
#' ambiguous strands ("*") are allowed to match with forward ("+") and reverse
#' ("-") strands.
#'
#' @param junc_df_all df that will contain info on junctions from all samples.
#' @param junc_df df containing the the info of junctions to be added.
#' @param i iteration of the loop. Used to give the count column a arbritrary,
#'   differentiating name for each sample. E.g. "count_1".
#'
#' @return df with the junctions from junc_df incoporated into junc_df_all.
.junc_merge <- function(junc_df_all, junc_df, i) {
    if (i == 1) {
        junc_df_all <- junc_df
    } else {

        # when merging allow for * strands to match with + or -
        junc_df_all <- junc_df_all %>%
            dplyr::full_join(junc_df,
                by = c("chr", "start", "end")
            ) %>%
            dplyr::mutate(
                strand = dplyr::case_when(
                    is.na(strand.y) ~ strand.x,
                    is.na(strand.x) ~ strand.y,
                    strand.x == strand.y ~ strand.x,
                    strand.x == "*" ~ strand.y,
                    strand.y == "*" ~ strand.x,
                    TRUE ~ NA_character_
                ),
                strand.x = NULL,
                strand.y = NULL
            )

        if (any(is.na(junc_df_all[["strand"]]))) {
            stop("No strands should be left as NA after processing")
        }
    }

    junc_df_all <- junc_df_all %>%
        dplyr::rename(!!stringr::str_c("count_", i) := count)

    return(junc_df_all)
}
