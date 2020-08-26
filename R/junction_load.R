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
#'   `junction_paths` with TRUE representing controls. Or, "fibroblasts"
#'   representing the samples of which GTEx tissue to use as controls. By
#'   default, will assume all samples are patients.
#' @param load_func function to load in junctions. By default, requires STAR
#'   formatted junctions (SJ.out). But this can be switched dependent on the
#'   format of the user's junction data. Function must take as input a junction
#'   path then return a data.frame with the columns "chr", "start", "end",
#'   "strand" and "count".
#' @param chrs chromosomes to keep. By default, no filter is applied.
#' @param coord_system One of "ensembl" (1-based) or "ucsc" (0-based) denoting
#'   the co-ordinate system corresponding to the user junctions from
#'   `junction_paths`. Only used when controls is set to "fibroblasts" to ensure
#'   GTEx data is harmonised to match the co-ordinate system of the user's
#'   junctions.
#'
#' @return
#' [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object containing junction data.
#'
#' @examples
#'
#' \dontrun{
#' # TO DO - figure out how to use system.file in examples
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
#' @family junction
#' @export
junction_load <- function(
    junction_paths,
    metadata = dplyr::tibble(samp_id = stringr::str_c("samp_", seq_along(junction_paths))),
    controls = rep(FALSE, length(junction_paths)),
    load_func = .STAR_load,
    chrs = NULL,
    coord_system = "ensembl") {
    # for R CMD Check
    chr <- NULL

    ##### Read in and merge junction data #####
    ##read as a list of data.table/data.frames
    junctions_list = lapply(junction_paths,load_func)
    ###add a SampleID column to each df in the list
    samp_id = stringr::str_c("samp_", seq_along(junction_paths)) ###HEADS UP not sure why metadata was a tibble??
    ###probably that can be removed
    junctions_list = purrr::map2(junctions_list, samp_id, ~cbind(.x, SampleID = .y))
    ###use the data.table rbindlist to make a long format
    junctions_all = data.table::rbindlist(junctions_list)
    if (!is.null(chrs)) {
      junctions_all <- .chr_filter(junctions_all, chrs) ###I think this should work fine
    }
    junctions_all <- .junction_merge(junctions_all) ###I think this should work fine


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
    # for R CMD Check
    chr <- `:=` <- uniq_map_read_count <- NULL #honestly no idea what going on here?
    ###using data.table to read in the files
    junctions <-
        data.table::fread(junction_path,
                          col.names = c(
                              "chr", "start", "end", "strand",
                              "intron_motif", "annotation", "uniq_map_read_count", "multi_map_read_count", "max_overhang"))
    ##data.table version of doing a case case_when
    ###we change type here so the as.character bit is required
    junctions[, strand := as.character(strand)][strand == "0", strand := "*"][strand == "1", strand := "+"][strand == "2", strand := "-"]
    ### .SD is subset data, .SDcols is either by name or number
    junctions = junctions[,.SD,.SDcols = c(1,2,3,4,7)]
    setnames(junctions,"uniq_map_read_count","count")


    return(junctions)
}

#' Merge two junction datasets together
#'
#' `.junction_merge` will merge a long dataset to wide format so will keep all rows from both datasets. It
#' will also ensure ambiguous strands ("*") are allowed to match with forward
#' ("+") and reverse ("-") strands.
#'
#' @param junctions_all data.table that will contain data on junctions from all in long
#'   samples.
#'
#' @return data.frame with the data `junctions_all` in a wide format
#'
#' @keywords internal
#' @noRd
.junction_merge <- function(junctions_all){
      ###wide to long not taking strand into account -- looks like this
    setkeyv(junctions_all, c("chr",  "start", "end","SampleID"))
    # chr start   end samp_1 samp_2 samp_3 samp_4 samp_5 samp_6 samp_7
    # 1: GL000195.1  2967 18392     NA     NA     NA     NA     NA      0     NA
    # 2: GL000195.1  2967 18609      1     NA     NA      1     NA     NA     NA
    # 3: GL000195.1 13868 18321     NA     NA     NA     NA     NA      0     NA
    junctions_counts = dcast(junctions_all, chr +  start + end ~ SampleID, value.var = c("count"))
    ###for every unique junciton take only the chr, start,end,strand
    strands = unique(junctions_all[,.SD,.SDcols = c(1,2,3,4)])
    ###now to a full join on the strands and counts tables
    setkeyv(junctions_counts, c("chr",  "start", "end"))
    setkeyv(strands, c("chr",  "start", "end"))
    with_strand = junctions_counts[strands,on = c("chr",  "start", "end")]
    ###now any start and end that had either '+' or '-' in one sample but "*" in another
    ### will be duplicated
    # chr start   end samp_1 samp_2 samp_3 samp_4 samp_5 samp_6 samp_7 strand
    # 1: GL000195.1  2967 18392     NA     NA     NA     NA     NA      0     NA      *
    #     2: GL000195.1  2967 18609      1     NA     NA      1     NA     NA     NA      *
    #     3: GL000195.1 13868 18321     NA     NA     NA     NA     NA      0     NA      *
    ###so we addd a count
    with_strand[,N_reps := .N, by = c("chr","start","end")]
    with_strand = with_strand[!(N_reps > 1 & strand == "*")]
    with_strand$N_reps = NULL

    return(with_strand)
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
