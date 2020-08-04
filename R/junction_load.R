#' Load junctions from patient and control RNA-seq data
#'
#' \code{junction_load} loads in raw patient and control junction data and
#' formats it into a
#' [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object. Control samples can be user-inputted or selected from GTEx data
#' publicly released through the recount2 project
#' (\url{https://jhubiostatistics.shinyapps.io/recount/}) and downloaded through
#' snaptron (\url{http://snaptron.cs.jhu.edu/}). By default,
#' \code{junction_load} expects the junction data to be in STAR aligned format
#' (SJ.out).
#'
#' @param junction_paths file path(s) to junction data.
#' @param metadata dataframe containing sample metadata with rows in the same
#'   order and corresponding to file path(s). Will be used as the \code{colData}
#'   of the
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   object.
#' @param controls either a logical vector of the same length as
#'   \code{junction_paths} with TRUEs labelling controls Or, "fibroblasts"
#'   representing the samples of which GTEx tissue to use as controls. If left
#'   unchanged, by default will assume all samples are patients.
#' @param load_func function to load in junctions. By default, requires STAR
#'   formatted junctions (SJ.out). But this can be switched dependent on the
#'   format of the user's junction data. Function must take as input a junction
#'   path then return a dataframe with the columns "chr", "start", "end",
#'   "strand" and "count".
#' @param chrs chromosomes to keep. If NULL, no filter is applied.
#' @param coord_system One of "ensembl" (1-based) or "ucsc" (0-based) denoting
#'   the co-ordinate system corresponding to the user junctions from
#'   \code{junction_paths}. Only used when controls is set to "fibroblasts".
#'   This is used ensure control data is harmonised to user's junctions when
#'   merging. The outputted junctions will always follow the user's co-ordinate
#'   system.
#'
#' @return
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   object containing junction data.
#'
#' @examples
#'
#' \dontrun{
#' junctions_example_1_path <-
#'     system.file("extdata", "junctions_example_1.txt",
#'         "dasper",
#'         mustWork = TRUE
#'     )
#'
#' junctions <-
#'     junction_load(
#'         junction_paths = c(junctions_example_1_path),
#'         metadata = dplyr::tibble(samp_id = c("example_1"))
#'     )
#'
#' junctions
#' }
#'
#' @export
junction_load <- function(junction_paths,
    metadata = dplyr::tibble(samp_id = stringr::str_c("samp_", seq_along(junction_paths))),
    controls = rep(FALSE, length(junction_paths)),
    load_func = .STAR_load,
    chrs = NULL,
    coord_system = "ensembl") {

    ##### Read in and merge junction data #####

    junctions_all <- NULL

    for (i in seq_along(junction_paths)) {
        print(stringr::str_c(Sys.time(), " - Loading junctions for sample ", i, "/", length(junction_paths), "..."))

        junctions <- load_func(junction_paths[i])

        if (!all(colnames(junctions) %in% c("chr", "start", "end", "strand", "count"))) {
            stop("load_func must return a dataframe with the columns 'chr', 'start', 'end', 'strand' and 'count'")
        }

        if (!is.null(chrs)) {
            junctions <- .chr_filter(junctions, chrs)
        }

        junctions_all <- .junction_merge(junctions_all, junctions)
    }

    ##### Add control data/identifier #####

    print(stringr::str_c(Sys.time(), " - Adding control junctions..."))

    if (any(controls %in% c("fibroblasts"))) {
        junctions_controls <- .junction_dl_controls(controls)

        junctions_controls <- .control_coord_convert(junctions_controls, coord_system)

        if (!is.null(chrs)) {
            junctions_controls <- .chr_filter(junctions_controls, chrs)
        }

        junctions_all <- .junction_merge(junctions_all, junctions_controls)

        metadata <- metadata %>%
            dplyr::mutate(case_control = c(rep("case", length(junction_paths)))) %>%
            dplyr::bind_rows(dplyr::tibble(case_control = rep("control", ncol(junctions_controls %>% dplyr::select(-chr, -start, -end, -strand)))))
    } else {
        if (!(is.logical(controls)) | length(controls) != length(junction_paths)) {
            stop("Controls argument must be a logical vector of same length as the number of junction paths or one of 'fibroblasts'")
        }

        metadata <- metadata %>%
            dplyr::mutate(case_control = ifelse(controls, "control", "case"))
    }

    ##### Tidy junction data #####

    print(stringr::str_c(Sys.time(), " - Tidying and storing junction data as a RangedSummarizedExperiment..."))

    # replace all missing count (NA) values with 0
    junctions_all[is.na(junctions_all)] <- 0

    # convert junctions into a RangedSummarizedExperiment
    raw_counts <- junctions_all %>%
        dplyr::select(-chr, -start, -end, -strand) %>%
        as.matrix()

    junction_coords <- junctions_all %>%
        dplyr::select(chr, start, end, strand) %>%
        GRanges()

    junctions <-
        SummarizedExperiment::SummarizedExperiment(
            assays = list(raw = raw_counts),
            rowRanges = junction_coords,
            colData = metadata
        )

    # sort ranges in natural order by chr, start, end
    junctions <- junctions %>%
        GenomeInfoDb::sortSeqlevels() %>%
        sort(ignore.strand = TRUE)

    print(stringr::str_c(Sys.time(), " - done!"))

    return(junctions)
}

#' Load raw junction data
#'
#' \code{.STAR_load} will load raw junction data that is outputted from STAR
#' (SJ.out) into R. This will format the junction data, retaining only chr,
#' start, end, strand and count (uniq_map_read_count) columns.
#'
#' @param junction_path path to the junction data.
#'
#' @return df detailing junction co-ordinates and counts.
#'
#' @keywords internal
#' @noRd
.STAR_load <- function(junction_path) {
    junctions <-
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

    return(junctions)
}

#' Merge two junction datasets together
#'
#' \code{.junction_merge} will merge two sets of junction data together. It uses a
#' full_join so will keep all rows from both datasets. It will also ensure
#' ambiguous strands ("*") are allowed to match with forward ("+") and reverse
#' ("-") strands.
#'
#' @param junctions_all df that will contain info on junctions from all samples.
#' @param junctions df containing the the info of junctions to be added.
#'
#' @return df with the junctions from junctions incoporated into junctions_all.
#'
#' @keywords internal
#' @noRd
.junction_merge <- function(junctions_all, junctions) {
    if (is.null(junctions_all)) {
        junctions_all <- junctions
    } else {

        # when merging allow for * strands to match with + or -
        junctions_all <- junctions_all %>%
            dplyr::full_join(junctions,
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

        if (any(is.na(junctions_all[["strand"]]))) {
            stop("No strands should be left as NA after processing")
        }
    }

    # convert count to count_X only if count exists as a column
    #
    if (any(colnames(junctions_all) == "count")) {

        # number of cols should equal the number of samples
        num_count_cols <- junctions_all %>%
            colnames() %>%
            stringr::str_detect("count") %>%
            sum()

        junctions_all <- junctions_all %>%
            dplyr::rename(!!stringr::str_c("count_", num_count_cols) := count)
    }

    return(junctions_all)
}

#' Download and merge control data from recount2
#'
#' \code{.junction_dl_controls} will download GTEx control junctions from
#' Dropbox using \code{\link[BiocFileCache]{BiocFileCache}}.
#'
#' @inheritParams junction_load
#'
#' @keywords internal
#' @noRd
.junction_dl_controls <- function(controls) {

    # generate details of controls to download
    controls_df <-
        dplyr::tibble(
            control = c("fibroblasts"),
            gtex_tissue = c("cells_transformed_fibroblasts"),
            dropbox_path = c("https://www.dropbox.com/s/bg54moy2qf2wnnm/GTEx_junctions_cells_transformed_fibroblasts.rda?dl=1")
        )

    controls_df <- controls_df %>%
        dplyr::filter(control == controls)

    print(stringr::str_c(
        Sys.time(),
        " - Downloading and importing ", controls_df[["control"]], " junction data..."
    ))

    file_path <- .file_cache(controls_df[["dropbox_path"]])

    load(file = file_path)

    return(GTEx_junctions_tidy)
}

#' Make sure that control junctions use the same coordinate system as the user's
#'
#' \code{.control_coord_convert} will convert control junctions to match users
#' junctions. If user's co-ordinate system set to "ensembl",  "chrM" will be
#' converted to "MT", the "chr" will be removed from chromosomes, and 1 will be
#' added to both "start" and "end". If user's co-ordinate system set to "ucsc",
#' 1 will be taken off from the start and end to convert from 1-based to 0-based
#' co-ordinates.
#'
#' @inheritParams junction_load
#'
#' @keywords internal
#' @noRd
.control_coord_convert <- function(junctions_controls, coord_system) {
    if (!any(coord_system %in% c("ensembl", "ucsc"))) {
        stop("coord_system must be one of 'ensembl' or 'ucsc'")
    }

    if (coord_system == "ensembl") {
        junctions_controls <- junctions_controls %>%
            dplyr::mutate(
                chr = chr %>% stringr::str_replace("chr", ""),
                chr = chr %>% stringr::str_replace("^M$", "MT")
            )
    } else if (coord_system == "ucsc") {
        junctions_controls <- junctions_controls %>%
            dplyr::mutate(
                start = start - 1,
                end = end - 1
            )
    }

    return(junctions_controls)
}
