#' Load junctions from patient and control RNA-seq data
#'
#' \code{junction_load} loads in raw patient and control junction data and formats
#' it into a \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' Control samples can be user-inputted or selected from GTEx data publicly
#' released through the recount2 project
#' (\url{https://jhubiostatistics.shinyapps.io/recount/}) and downloaded
#' through snaptron (\url{http://snaptron.cs.jhu.edu/}). By default, \code{junction_load}
#' expects the junction data to be in STAR aligned format (SJ.out).
#'
#' @param junction_paths file path(s) to junction data.
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
#' @param chrs chrs chromosomes to keep. If NULL, no filter is applied.
#'
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   containing junction data.
#'
#' @examples
#'
#' \dontrun{
#' example_junctions_1_path <-
#'     system.file("extdata", "example_junctions_1.txt",
#'         "dasper",
#'         mustWork = TRUE
#'     )
#'
#' junctions <-
#'     junction_load(
#'         junction_paths = c(example_junctions_1_path),
#'         metadata = dplyr::tibble(samp_id = c("example_1"))
#'     )
#'
#' junctions
#' }
#'
#' @export
junction_load <- function(junction_paths,
    metadata = NULL,
    controls = rep(FALSE, length(junction_paths)),
    load_func = .load_STAR,
    chrs = NULL) {

    ##### Read in and merge junction data #####

    junction_df_all <- NULL

    for (i in seq_along(junction_paths)) {
        print(stringr::str_c(Sys.time(), " - loading junctions for sample ", i, "/", length(junction_paths), "..."))

        junction_df <- load_func(junction_paths[i])

        if (!all(colnames(junction_df) %in% c("chr", "start", "end", "strand", "count"))) {
            stop("load_func must return a dataframe with the columns 'chr', 'start', 'end', 'strand' and 'count'")
        }

        if (!is.null(chrs)) {
            junction_df <- .chr_filter(junction_df, chrs)
        }

        junction_df_all <- .junction_merge(junction_df_all, junction_df)
    }

    ##### Add control data/identifier #####

    if (any(controls %in% c("fibroblasts", "whole_blood"))) {



    } else {
        if (!(is.logical(controls)) | length(controls) != length(junction_paths)) {
            stop("Controls argument must be a logical vector of same length as the number of junction paths")
        }

        metadata <- metadata %>%
            dplyr::mutate(case_control = ifelse(controls, "control", "case"))
    }

    ##### Tidy junction data #####

    # replace all missing count (NA) values with 0
    junction_df_all[is.na(junction_df_all)] <- 0

    # convert junctions into a RangedSummarizedExperiment
    raw_counts <- junction_df_all %>%
        dplyr::select(-chr, -start, -end, -strand) %>%
        as.matrix()

    junction_coords <- junction_df_all %>%
        dplyr::select(chr, start, end, strand) %>%
        GRanges()

    junctions <-
        SummarizedExperiment::SummarizedExperiment(
            assays = list(raw_counts),
            rowRanges = junction_coords,
            colData = metadata
        )

    print(stringr::str_c(Sys.time(), " - done!"))

    return(junctions)
}

#' Load raw junction data
#'
#' \code{.load_STAR} will load raw junction data that is outputted from STAR
#' (SJ.out) into R. This will format the junction data, retaining only chr,
#' start, end, strand and count (uniq_map_read_count) columns.
#'
#' @param junction_path path to the junction data.
#'
#' @return df detailing junction co-ordinates and counts.
#'
#' @keywords internal
#' @noRd
.load_STAR <- function(junction_path) {
    junction_df <-
        readr::read_delim(junction_path,
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

    return(junction_df)
}



#' Merge two junction datasets together
#'
#' \code{.junction_merge} will merge two sets of junction data together. It uses a
#' full_join so will keep all rows from both datasets. It will also ensure
#' ambiguous strands ("*") are allowed to match with forward ("+") and reverse
#' ("-") strands.
#'
#' @param junction_df_all df that will contain info on junctions from all samples.
#' @param junction_df df containing the the info of junctions to be added.
#' @param i iteration of the loop. Used to give the count column a arbritrary,
#'   differentiating name for each sample. E.g. "count_1".
#'
#' @return df with the junctions from junction_df incoporated into junction_df_all.
#'
#' @keywords internal
#' @noRd
.junction_merge <- function(junction_df_all, junction_df) {
    if (is.null(junction_df_all)) {
        junction_df_all <- junction_df
    } else {

        # when merging allow for * strands to match with + or -
        junction_df_all <- junction_df_all %>%
            dplyr::full_join(junction_df,
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

        if (any(is.na(junction_df_all[["strand"]]))) {
            stop("No strands should be left as NA after processing")
        }
    }

    # number of cols should equal the number of samples
    # use this value to annotate the co
    num_count_cols <- junction_df_all %>%
        colnames() %>%
        stringr::str_detect("count") %>%
        sum()

    junction_df_all <- junction_df_all %>%
        dplyr::rename(!!stringr::str_c("count_", num_count_cols) := count)

    return(junction_df_all)
}

#' Download and merge control data from recount2
#'
#' \code{.junction_add_controls}
#'
#' @keywords internal
#' @noRd
.junction_add_controls <- function(controls, junction_df_all, dest_dir = tempdir()) {

    # generate details of controls to download
    controls_df <-
        dplyr::tibble(
            control = c("fibroblasts"),
            gtex_tissue = c("cells_transformed_fibroblasts"),
            dropbox_path = c("https://www.dropbox.com/s/6w3nrbzt3nknkmh/GTEx_junctions_cells_transformed_fibroblasts.rda?dl=1"),
            file_name = dropbox_path %>%
                stringr::str_replace(".*/", "") %>%
                stringr::str_replace("\\?.*", "")
        )

    controls_df <- controls_df %>%
        dplyr::filter(control == controls)

    # create a temporary destination folder and download control junction data
    dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

    file_path <- stringr::str_c(dest_dir, "/", controls_df[["file_name"]])

    print(stringr::str_c(Sys.time(), " - Downloading junction data from GTEx controls: ", controls_df[["control"]], "..."))

    utils::download.file(
        url = controls_df[["dropbox_path"]],
        destfile = file_path,
        mode = "wb"
    )

    load(file = file_path)

    # merge control data with cases
    junction_df_all <- junction_df_all %>%
        .junction_merge(GTEx_junctions_tidy)

    return(junction_df_all)
}
