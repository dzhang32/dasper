#' Load junctions from RNA-sequencing data
#'
#' `junction_load` loads in raw patient and control junction data and formats it
#' into a
#' [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object. Control samples can be the user's in-house samples or selected from
#' GTEx v6 data publicly released through the
#' \href{https://jhubiostatistics.shinyapps.io/recount/}{recount2} and
#' downloaded through \href{http://snaptron.cs.jhu.edu/}{snaptron}. By default,
#' `junction_load` expects the junction data to be in STAR aligned format
#' (SJ.out) but this can be modified via the argument `load_func`.
#'
#' @param junction_paths path(s) to junction data.
#' @param metadata data.frame containing sample metadata with rows in the same
#'   order as `junction_paths`.
#' @param controls either a logical vector of the same length as
#'   `junction_paths` with TRUE representing controls. Or, one of "fibroblasts",
#'   "lymphocytes", "skeletal_muscle", "whole_blood" representing the samples of
#'   which GTEx tissue to use as controls. By default, will assume all samples
#'   are patients.
#' @param load_func function to load in junctions. By default, requires STAR
#'   formatted junctions (SJ.out). But this can be switched dependent on the
#'   format of the user's junction data. Function must take as input a junction
#'   path then return a data.frame with the columns "chr", "start", "end",
#'   "strand" and "count".
#' @param chrs chromosomes to keep. By default, no filter is applied.
#' @param coord_system 1 (1-based) or 0 (0-based) denoting the co-ordinate
#'   system corresponding to the user junctions from `junction_paths`. Only used
#'   when controls is set to "fibroblasts" to ensure GTEx data is harmonised to
#'   match the co-ordinate system of the user's junctions.
#'
#' @return
#' [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object containing junction data.
#'
#' @examples
#'
#' junctions_example_1_path <-
#'     system.file("extdata",
#'         "junctions_example_1.txt",
#'         package = "dasper",
#'         mustWork = TRUE
#'     )
#' junctions_example_2_path <-
#'     system.file("extdata",
#'         "junctions_example_2.txt",
#'         package = "dasper",
#'         mustWork = TRUE
#'     )
#'
#' junctions <-
#'     junction_load(
#'         junction_paths = c(junctions_example_1_path, junctions_example_2_path)
#'     )
#'
#' junctions
#' @export
junction_load <- function(
    junction_paths,
    metadata = dplyr::tibble(samp_id = stringr::str_c("samp_", seq_along(junction_paths))),
    controls = rep(FALSE, length(junction_paths)),
    load_func = .STAR_load,
    chrs = NULL,
    coord_system = 1) {
    # for R CMD Check
    chr <- NULL

    ##### Read in and merge junction data #####

    junctions_all <- NULL

    for (i in seq_along(junction_paths)) {
        print(stringr::str_c(Sys.time(), " - Loading junctions for sample ", i, "/", length(junction_paths), "..."))

        junctions <- load_func(junction_paths[i])

        if (!all(colnames(junctions) %in% c("chr", "start", "end", "strand", "count"))) {
            stop("load_func must return a data.frame with the columns 'chr', 'start', 'end', 'strand' and 'count'")
        }

        if (!is.null(chrs)) {
            junctions <- .chr_filter(junctions, chrs)
        }

        junctions_all <- .junction_merge(junctions_all, junctions)
    }

    ##### Add control data/identifier #####

    print(stringr::str_c(Sys.time(), " - Adding control junctions..."))

    if (any(controls %in% c("fibroblasts", "lymphocytes", "skeletal_muscle", "whole_blood"))) {
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

#' Load STAR-formatted junctions
#'
#' `.STAR_load` will load raw junction data that is outputted from STAR
#' (SJ.out) and format it, retaining only chr, start, end, strand and count
#' (uniq_map_read_count) columns.
#'
#' @param junction_path path to the junction data.
#'
#' @return data.frame detailing junction co-ordinates and counts.
#'
#' @keywords internal
#' @noRd
.STAR_load <- function(junction_path) {

    # using data.table to read in the files
    junctions <- data.table::fread(junction_path,
        col.names = c(
            "chr", "start", "end", "strand",
            "intron_motif", "annotation", "uniq_map_read_count", "multi_map_read_count", "max_overhang"
        )
    )

    # data.table version of doing a case case_when
    # we change type here so the as.character bit is required
    junctions[, strand := as.character(strand)][strand == "0", strand := "*"][strand == "1", strand := "+"][strand == "2", strand := "-"]

    # .SD is subset data, .SDcols is either by name or number
    junctions <- junctions[, .SD, .SDcols = c(1, 2, 3, 4, 7)]
    data.table::setnames(junctions, "uniq_map_read_count", "count")

    junctions <- junctions %>% as.data.frame()

    return(junctions)
}

#' Merge two junction datasets together
#'
#' `.junction_merge` will merge two sets of junction data together. It uses a
#' [full_join][dplyr::mutate-joins] so will keep all rows from both datasets. It
#' will also ensure ambiguous strands ("*") are allowed to match with forward
#' ("+") and reverse ("-") strands.
#'
#' @param junctions_all data.frame that will contain data on junctions from all
#'   samples.
#' @param junctions data.frame containing the data of junctions from one sample.
#'
#' @return data.frame with the data `junctions` incoporated into
#'   `junctions_all`.
#'
#' @keywords internal
#' @noRd
.junction_merge <- function(junctions_all, junctions) {

    # for R CMD Check
    `:=` <- NULL

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
#' @return data.frame containing GTEx junctions.
#'
#' @keywords internal
#' @noRd
.junction_dl_controls <- function(controls) {

    # for R CMD check
    control <- GTEx_junctions_tidy <- NULL

    # generate details of controls to download
    controls_df <-
        dplyr::tibble(
            control = c("fibroblasts", "lymphocytes", "skeletal_muscle", "whole_blood"),
            gtex_tissue = c(
                "cells_transformed_fibroblasts",
                "lymphocytes",
                "skeletal_muscle", "
                            whole_blood"
            ),
            dropbox_path = c(
                "https://www.dropbox.com/s/kh9fvknltai89fo/GTEx_v6_junctions_fibroblast.rda?dl=1",
                "https://www.dropbox.com/s/trtexmcuw9m8whj/GTEx_v6_junctions_lymphocytes.rda?dl=1",
                "https://www.dropbox.com/s/yo8hrqs75l4q03m/GTEx_v6_junctions_skeletal_muscle.rda?dl=1",
                "https://www.dropbox.com/s/hfd1go6ynztsjdy/GTEx_v6_junctions_whole_blood.rda?dl=1"
            )
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

#' Harmonise the control junction co-ordinates
#'
#' `.control_coord_convert` will convert control junctions to match users
#' junctions. If `coord_system` set to "ensembl", "chrM" will be converted to
#' "MT", the "chr" will be removed from chromosomes. If `coord_system` set to
#' "ucsc", 1 will be taken off from the start and end to convert from 1-based to
#' 0-based co-ordinates.
#'
#' @inheritParams junction_load
#'
#' @return data.frame containing control junctions with harmonised co-ordinates.
#'
#' @keywords internal
#' @noRd
.control_coord_convert <- function(junctions_controls, coord_system) {

    # for R CMD check
    chr <- NULL

    if (!any(coord_system %in% c(1, 0))) {
        stop("coord_system must be either 1 or 0")
    }

    if (coord_system == 0) {
        junctions_controls <- junctions_controls %>%
            dplyr::mutate(
                start = start - 1,
                end = end - 1
            )
    }

    return(junctions_controls)
}
