#' Visualise RNA-seq data in a the form of a sashimi plot
#'
#' \code{plot_sashimi} plots the splicing events and coverage across specific
#' genes/transcripts/regions of interest. Unlike traditional sashimi plots,
#' coverage and junction tracks are separated, which enables user's to choose
#' whether they would like to plot only the junctions.
#'
#' @inheritParams junction_annot
#' @inheritParams coverage_process
#'
#' @param gene_tx_id character scalar with the id of the gene. Currently, this
#'   must be in the form of an Ensembl gene or transcript id, which has a
#'   matching entry in `ref`.
#' @param case_id list containing 1 element. The contents of this element must
#'   be a character vector specifying sample ids that are to be plotted. The
#'   name of this element must correspond to the column containing sample ids in
#'   the junction `SummarizedExperiment::mcols()`. By default, all cases will be
#'   plotted.
#' @param sum_func function that will be used to aggregate the junction counts
#'   and coverage for controls. By default, `mean` will be used.
#' @param region a [GenomicRanges][GenomicRanges::GRanges-class] of length 1
#'   that is used to filter the exons/junctions plotted. Only those that overlap
#'   this region are plotted.
#' @param annot_colour character vector length 7, representing the colours of
#'   junction types.
#' @param digits used in `round`, specifying the number of digits to round the
#'   junction counts to for visualisation purposes.
#' @param count_label logical value specifying whether to add label the count of
#'   each junction.
#' @param load_func function used to load coverage.
#' @param binwidth the number of bases to aggregate coverage across using
#'   `sum_func` when plotting. .
#'
#' @return `ggplot` displaying the splicing (and coverage) surrounding the
#'   transcript/region of interest.
#'
#' @examples
#'
#' # use GenomicState to load txdb (GENCODE v31)
#' ref <- GenomicState::GenomicStateHub(
#'     version = "31",
#'     genome = "hg38",
#'     filetype = "TxDb"
#' )[[1]]
#'
#' junctions_processed <- junction_process(
#'     junctions_example,
#'     ref,
#'     types = c("ambig_gene", "unannotated")
#' )
#'
#' sashimi_plot <- plot_sashimi(
#'     junctions = junction_filter(junctions_processed),
#'     ref = ref,
#'     gene_tx_id = "ENSG00000142156.14",
#'     sum_func = NULL
#' )
#' @export
plot_sashimi <- function(junctions,
    ref,
    gene_tx_id,
    case_id = NULL,
    sum_func = mean,
    region = NULL,
    annot_colour = c(
        ggpubr::get_palette("jco", 1),
        ggpubr::get_palette("npg", 7)[c(1, 3, 2, 5, 6)],
        ggpubr::get_palette("jco", 6)[c(3)]
    ),
    digits = 2,
    count_label = TRUE,
    coverage_paths_case = NULL,
    coverage_paths_control = NULL,
    coverage_chr_control = NULL,
    load_func = .coverage_load,
    binwidth = 100) {

    ##### Load reference annotation #####

    ref <- .ref_load(ref)

    ##### Obtain the exons and junctions to plot #####

    gene_tx_list <- .gene_tx_type_get(gene_tx_id)

    exons_to_plot <- .exons_to_plot_get(ref, gene_tx_list, region)

    junctions_to_plot <- .junctions_to_plot_get(junctions, gene_tx_list, region)

    ##### Obtain co-ordinates to plot #####

    gene_tx_to_plot <- GenomicFeatures::genes(ref, filter = gene_tx_list)

    coords_to_plot <- .coords_to_plot_get(gene_tx_to_plot, exons_to_plot, junctions_to_plot)

    ##### Plot gene (exons and gene body line) #####

    gene_track_plot <- .plot_gene_track(coords_to_plot, exons_to_plot)

    ##### Plot junctions #####

    sashimi_plots <- .plot_sashimi_junctions(
        junctions_to_plot,
        gene_track_plot,
        case_id,
        sum_func,
        digits,
        assay_name = "norm",
        annot_colour,
        count_label
    )

    ##### Plot coverage #####

    if (!is.null(coverage_paths_case)) {
        coverage_to_plot <- .coverage_to_plot_get(coords_to_plot,
            coverage_paths_case,
            coverage_paths_control,
            coverage_chr_control = coverage_chr_control,
            load_func = .coverage_load,
            sum_func = sum_func
        )

        coverage_plot <- .plot_coverage(
            coverage_to_plot,
            coords_to_plot,
            binwidth
        )

        sashimi_plots <- .merge_coverage_sashimi(coverage_plot, sashimi_plots)
    }

    ##### Arrange plots #####

    sashimi_plots <- ggpubr::ggarrange(
        plotlist = sashimi_plots,
        ncol = 1,
        nrow = length(sashimi_plots),
        common.legend = TRUE,
        align = "v",
        legend = "top"
    )

    ##### Add annotation #####

    sashimi_plots <- .plot_annotation(sashimi_plots, gene_tx_id, coords_to_plot)

    return(sashimi_plots)
}

#' Tidy user-inputted gene or transcript id
#'
#' `.gene_tx_type_get` will recognize whether the user has inputted a gene or
#' transcript ID. Then derive the column which the gene/transcript should be
#' matched against in `junctions`.
#'
#' @param gene_tx_id gene or transcript id (currently MUST be in a "ENSG" or
#'   "ENST" format).
#'
#' @return list with gene/transcript id, the name of which corresponds to the
#'   `SummarizedExperiment::rowRanges` column to filter in `junctions`.
#'
#' @keywords internal
#' @noRd
.gene_tx_type_get <- function(gene_tx_id) {
    if (is.null(gene_tx_id) | length(gene_tx_id) != 1) {
        stop("gene_tx_id must be set and be of length 1")
    } else if (stringr::str_detect(gene_tx_id, "ENSG")) {
        gene_tx_type <- "gene_id"
    } else if (stringr::str_detect(gene_tx_id, "ENST")) {
        gene_tx_type <- "tx_name"
    } else {
        stop("gene_tx_id does not include an ENST or ENSG prefix")
    }

    # create named list of gene/tx id
    # for filtering txdb
    gene_tx_list <- list(gene_tx_id)
    names(gene_tx_list) <- gene_tx_type

    return(gene_tx_list)
}

#' Obtain exons to be plotted
#'
#' `.exons_to_plot_get` will obtain the exons from the reference annotation base
#' on the user-inputted gene or transcript of interest. Then will filter for
#' only exons overlapping `region`. Will use `GenomicRanges::disjoin` to
#' collapse together overlapping exons.
#'
#' @inheritParams plot_sashimi
#'
#' @param gene_tx_list list containing gene/transcript id returned from
#'   `.gene_tx_type_get`.
#'
#' @return [GenomicRanges][GenomicRanges::GRanges-class] object containing exons
#'   to be plotted.
#'
#' @keywords internal
#' @noRd
.exons_to_plot_get <- function(ref,
    gene_tx_list,
    region) {

    # filter for exons of gene/tx of interest
    exons_to_plot <- GenomicFeatures::exons(ref, filter = gene_tx_list)

    # if exons overlap (e.g. for gene inputs), disjoin for plotting
    exons_to_plot <- exons_to_plot %>% GenomicRanges::disjoin()

    # keep only exons that overlap your region of interest
    if (!is.null(region)) {
        if (!is(region, "GenomicRanges") | length(region) != 1) {
            stop("region must be a GenomicRanges object of length 1")
        }

        exon_region_hits <- findOverlaps(query = exons_to_plot, subject = region)
        exons_to_plot <- exons_to_plot[S4Vectors::queryHits(exon_region_hits)]
    }

    if (length(exons_to_plot) == 0) {
        stop("No exons found to plot")
    }

    return(exons_to_plot)
}

#' Obtain junctions to be plotted
#'
#' `.junctions_to_plot_get` will obtain the junctions that are connected with
#' the gene or transcript via `junction_annot` and fall within `region`.
#'
#' @inheritParams plot_sashimi
#'
#' @return
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   object containing junctions to be plotted.
#'
#' @keywords internal
#' @noRd
.junctions_to_plot_get <- function(junctions, gene_tx_list, region) {
    gene_tx <- gene_tx_list %>% unlist()

    # check the columns used are in a CharacterList format
    col_type_chr_list <-
        methods::is(GenomicRanges::mcols(junctions)[[stringr::str_c(names(gene_tx), "_start")]], "CharacterList") &
            methods::is(GenomicRanges::mcols(junctions)[[stringr::str_c(names(gene_tx), "_end")]], "CharacterList")

    if (!col_type_chr_list) {
        stop(stringr::str_c(
            "Columns storing the gene/tx information are not CharacterList objects",
            " - was this SE generated using dasper::junction_annot()?"
        ))
    }

    # currently the any() used here may be too liberal, especially for overlapping genes
    # but may be okay, since juncs need to precisely match the exon boundary
    junctions_indexes <-
        which(any(GenomicRanges::mcols(junctions)[[stringr::str_c(names(gene_tx), "_start")]] == gene_tx) |
            any(GenomicRanges::mcols(junctions)[[stringr::str_c(names(gene_tx), "_end")]] == gene_tx))

    junctions_to_plot <- junctions[junctions_indexes, ]

    # keep only junctions that overlap your region of interest
    if (!is.null(region)) {
        junctions_region_hits <- findOverlaps(query = junctions_to_plot, subject = region)
        junctions_to_plot <- junctions_to_plot[S4Vectors::queryHits(junctions_region_hits)]
    }

    if (length(junctions_to_plot) == 0) {
        stop("No junctions found to plot")
    }

    return(junctions_to_plot)
}

#' Obtains the co-ordinates used for plotting
#'
#' `.coords_to_plot_get` will obtain the junctions that are connected with the
#' gene or transcript via `junction_annot` and fall within `region`.
#'
#' @inheritParams plot_sashimi
#'
#' @param exons_to_plot [GenomicRanges][GenomicRanges::GRanges-class] object
#'   containing exons to be plotted.
#' @param junctions_to_plot
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'    object containing junctions to be plotted.
#' @param ext_factor the factor by which to extend the x-axes limits of the
#'   plot.
#'
#' @return list containing the coordinates to be used for plotting.
#'
#' @keywords internal
#' @noRd
.coords_to_plot_get <- function(gene_tx_to_plot, exons_to_plot, junctions_to_plot, ext_factor = 30) {
    coords_to_plot <-
        list(
            chr = gene_tx_to_plot %>% seqnames() %>% as.character() %>% unique(),
            strand = gene_tx_to_plot %>% strand() %>% as.character() %>% unique(),
            gene_start = gene_tx_to_plot %>% start(),
            gene_end = gene_tx_to_plot %>% end(),
            min_exon_start = c(start(exons_to_plot), start(junctions_to_plot)) %>% min(),
            max_exon_end = c(end(exons_to_plot), end(junctions_to_plot)) %>% max()
        )

    # add a gap between end of exon and end of plot for visualisation
    coords_to_plot[["range_exon_start_end"]] <- coords_to_plot[["max_exon_end"]] - coords_to_plot[["min_exon_start"]]
    ext_num <- round(coords_to_plot[["range_exon_start_end"]] / ext_factor)
    coords_to_plot[["x_min"]] <- coords_to_plot[["min_exon_start"]] - ext_num
    coords_to_plot[["x_max"]] <- coords_to_plot[["max_exon_end"]] + ext_num

    # set coords of the line segment marking gene body
    # extend this if end of gene/tx falls outside of x_min/x_max
    coords_to_plot[["segment_start"]] <- ifelse(coords_to_plot[["gene_start"]] < coords_to_plot[["x_min"]],
        coords_to_plot[["x_min"]],
        coords_to_plot[["min_exon_start"]]
    )

    coords_to_plot[["segment_end"]] <- ifelse(coords_to_plot[["gene_end"]] > coords_to_plot[["x_max"]],
        coords_to_plot[["x_max"]],
        coords_to_plot[["max_exon_end"]]
    )

    # convert numeric values to integers
    coords_to_plot <- coords_to_plot %>%
        lapply(FUN = function(x) if (is.numeric(x)) as.integer(x) else x)

    return(coords_to_plot)
}

#' Plot the exons and gene/transcript of interest
#'
#' `.plot_gene_track` will plot the exons and gene body of the inputted gene as a `ggplot`.
#'
#' @inheritParams plot_sashimi
#'
#' @param coords_to_plot list containing the coordinates to be used for plotting.
#' @param exons_to_plot [GenomicRanges][GenomicRanges::GRanges-class] object
#'   containing exons to be plotted.
#'
#' @return A `ggplot` object with the gene track.
#'
#' @keywords internal
#' @noRd
.plot_gene_track <- function(coords_to_plot, exons_to_plot) {

    # convert to df for ggplot
    exons_to_plot_df <- exons_to_plot %>% as.data.frame()

    # plot gene line
    gene_track <- ggplot2::ggplot() +
        ggplot2::geom_segment(ggplot2::aes(
            x = coords_to_plot[["segment_start"]],
            xend = coords_to_plot[["segment_end"]],
            y = 0, yend = 0
        ),
        size = 2
        )

    # plot exons
    gene_track <- gene_track +
        ggplot2::geom_rect(
            data = exons_to_plot_df,
            ggplot2::aes(
                xmin = start, xmax = end,
                ymin = -0.25, ymax = 0.25
            ),
            colour = "black",
            fill = ggpubr::get_palette("Greys", 10)[4]
        )

    # plot strand arrow
    gene_track <- gene_track +
        ggplot2::geom_segment(ggplot2::aes(
            x = ifelse(coords_to_plot[["strand"]] == "+",
                coords_to_plot[["segment_start"]],
                coords_to_plot[["segment_end"]]
            ),
            xend = ifelse(coords_to_plot[["strand"]] == "+",
                coords_to_plot[["segment_start"]] + coords_to_plot[["range_exon_start_end"]] / 30,
                coords_to_plot[["segment_end"]] - coords_to_plot[["range_exon_start_end"]] / 30
            ),
            y = 0.65, yend = 0.65
        ),
        size = 0.75,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.3, units = "cm"))
        )

    # add scale/theme aesthetic tweaks
    gene_track <- gene_track +
        ggplot2::scale_y_continuous(limits = c(-1, 1)) +
        ggplot2::scale_x_continuous(
            name = stringr::str_c("Chromosome ", coords_to_plot[["chr"]]),
        ) +
        ggplot2::coord_cartesian(
            xlim = c(
                coords_to_plot[["x_min"]],
                coords_to_plot[["x_max"]]
            )
        ) +
        ggpubr::theme_pubclean(flip = TRUE) +
        ggplot2::theme(
            axis.line.y = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank()
        )

    return(gene_track)
}

#' Plot the junction track for cases and controls
#'
#' `.plot_sashimi_junctions` will plot the junctions overlayed on the
#' `gene_track_plot` for cases and controls. Internally, this uses
#' `.junctions_counts_type_get` and `.junctions_points_get` and
#' `.plot_junctions`.
#'
#' @inheritParams plot_sashimi
#'
#' @param junctions_to_plot
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'    object containing junctions to be plotted.
#' @param gene_track_plot `ggplot` object displaying the exons and gene body
#'   returned by `.plot_gene_track`.
#' @param assay_name character scalar specifying the
#'   `SummarizedExperiment::assay()` from which to obtain junction counts.
#'
#' @return A `ggplot` object with the exons and junctions for cases and
#'   controls.
#'
#' @keywords internal
#' @noRd
.plot_sashimi_junctions <- function(junctions_to_plot,
    gene_track_plot,
    case_id,
    sum_func,
    digits,
    assay_name,
    annot_colour,
    count_label) {

    # format junctions into df with count/type details
    junctions_to_plot <- .junctions_counts_type_get(
        junctions_to_plot = junctions_to_plot,
        case_id = case_id,
        sum_func = sum_func,
        digits = digits,
        assay_name = assay_name
    )

    # obtain points for curves of junctions
    junctions_to_plot <- .junctions_points_get(
        junctions_counts = junctions_to_plot,
        ncp = 25
    )

    sashimi_plots <- .plot_junctions(
        junctions_to_plot = junctions_to_plot,
        gene_track_plot = gene_track_plot,
        annot_colour = annot_colour,
        count_label = count_label
    )

    return(sashimi_plots)
}

#' Obtain the junction counts
#'
#' `.junctions_counts_type_get` will obtain the counts to plot for the
#' cases/controls. For controls, if there is more than one sample, it will
#' summarise the counts via `sum_func`.
#'
#' @inheritParams plot_sashimi
#'
#' @param junctions_to_plot
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'    object containing junctions to be plotted.
#'
#' @return A `data.frame` will the junctions to be plotted and their associated
#'   counts for cases/controls.
#'
#' @keywords internal
#' @noRd
.junctions_counts_type_get <- function(junctions_to_plot,
    case_id = list(samp_id = "samp_1"),
    sum_func = mean,
    digits = 2,
    assay_name = "norm") {

    # for R CMD Check
    index <- type <- . <- NULL

    # check that the assay exists
    if (!(assay_name %in% SummarizedExperiment::assayNames(junctions_to_plot))) {
        stop(stringr::str_c(assay_name, " not found in junctions object"))
    }

    junctions_counts <- junctions_to_plot %>%
        GenomicRanges::ranges() %>%
        as.data.frame() %>%
        dplyr::mutate(
            index = dplyr::row_number(),
            type = mcols(junctions_to_plot)[["type"]]
        )

    # make sure we only have the required columns
    # specfically, avoid error induced by rownames of a GRanges
    junctions_counts <- junctions_counts %>%
        dplyr::select(start, end, width, index, type)

    if (is.null(case_id)) {
        samp_ids <- stringr::str_c("samp_", c(seq_len(dim(junctions_to_plot)[2])))
    } else {
        samp_id_col <- names(case_id)
        samp_ids <- case_id[[samp_id_col]]
    }

    # retrieve counts for the samples of interest
    for (i in seq_along(samp_ids)) {
        if (is.null(case_id)) {
            samp_index <- i
        } else {
            samp_index <- which(colData(junctions_to_plot)[[samp_id_col]] == samp_ids[i])
        }

        junctions_counts[[samp_ids[i]]] <- assays(junctions_to_plot)[[assay_name]][, samp_index] %>%
            round(digits = digits)
    }

    # aggregate and add counts for control samples
    if (!is.null(sum_func)) {
        which_control <- which(colData(junctions_to_plot)[["case_control"]] == "control")

        junctions_counts[["control"]] <-
            assays(junctions_to_plot)[[assay_name]][, which_control] %>%
            apply(MARGIN = 1, FUN = sum_func) %>%
            round(digits = digits)
    }

    # filter for out junctions that are not expressed (> 0 counts) in any sample
    junctions_counts <- junctions_counts %>%
        dplyr::select(-start, -end, -width, -index, -type) %>%
        apply(MARGIN = 1, FUN = function(x) !all(x == 0)) %>%
        dplyr::filter(junctions_counts, .)

    return(junctions_counts)
}

#' Obtain the junction points
#'
#' `.junctions_points_get` will obtain the points of the arc used to plot each
#' junction. This is based off of the `grid:::calcControlPoints()` function.
#' Will also mark the midpoint of each junction, used for the labelled of
#' junction counts.
#'
#' @param junctions_counts `data.frame` containing junction counts returned by
#'   `.junctions_counts_type_get`.
#'
#' @return A `data.frame` with junction x-coords calculated.
#'
#' @keywords internal
#' @noRd
.junctions_points_get <- function(junctions_counts, ncp = 25) {
    # For R CMD Check
    y <- index <- x <- . <- type <- mid_point <- samp_id <- NULL

    # calculate the points to plot the curve for each junction
    # without this (using ggplot2::geom_curve), there's no way of knowing the midpoint of y to add a label
    junctions_points <- grid:::calcControlPoints(
        x1 = junctions_counts[["start"]], x2 = junctions_counts[["end"]],
        y1 = 0, y2 = 0,
        angle = 90,
        curvature = -0.5,
        ncp = ncp
    ) # how many ctrl points per junc?

    junctions_points <-
        dplyr::tibble(
            x = junctions_points[["x"]],
            y = junctions_points[["y"]],
            index = junctions_counts[["index"]] %>%
                rep(times = ncp) %>% # repeat these indexes the same number of ctrl points
                sort()
        )

    # add the intial start/end positions, since these are not return from the ctrl points
    junctions_points <-
        dplyr::tibble(
            x = c(junctions_counts[["start"]], junctions_counts[["end"]]),
            y = 0, # start from middle of exons
            index = rep(junctions_counts[["index"]], 2)
        ) %>%
        dplyr::bind_rows(junctions_points)

    # normalise y so values should always sit between -1 and 1
    # change, even junctions y values to negative - to be plotted on bottom
    junctions_points <- junctions_points %>%
        dplyr::mutate(
            y = y / max(y),
            y = ifelse(index %% 2 == 0, -y, y)
        )

    # mark midpoints of junction curves to add label
    junctions_points <- junctions_points %>%
        dplyr::group_by(index) %>%
        dplyr::mutate(mid_point = .mark_mid(x)) %>%
        dplyr::ungroup()

    # add junction categories to colour by
    # add counts and tidy from wide into a long format for plotting
    junctions_points <- junctions_counts %>%
        dplyr::select(-start, -end, -width) %>%
        dplyr::left_join(junctions_points, ., by = "index") %>%
        dplyr::mutate(type = type %>% factor(
            levels = c(
                "annotated",
                "novel_acceptor",
                "novel_donor",
                "novel_combo",
                "novel_exon_skip",
                "ambig_gene",
                "unannotated"
            )
        )) %>%
        tidyr::gather(
            key = "samp_id", value = "count",
            -index, -type, -x, -y, -mid_point
        )

    junctions_points <- junctions_points %>%
        dplyr::mutate(linetype = ifelse(count == 0, 2, 1)) %>%
        dplyr::arrange(samp_id, index, x)

    return(junctions_points)
}

#' Plot the junction counts for a particular sample
#'
#' `.plot_junctions` will plot the junctions for each sample and store outputted
#' `ggplot`s into a list.
#'
#' @param junctions_to_plot
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'    object containing junctions to be plotted.
#' @param gene_track_plot `ggplot` object displaying the exons and gene body
#'   returned by `.plot_gene_track`.
#'
#' @return A `list` containing 1 `ggplot` object per sample.
#'
#' @keywords internal
#' @noRd
.plot_junctions <- function(junctions_to_plot = junctions_to_plot,
    gene_track_plot = gene_track_plot,
    annot_colour = annot_colour,
    count_label = count_label) {
    # for R CMD Check
    samp_id <- x <- y <- index <- type <- mid_point <- NULL
    samp_ids <- unique(junctions_to_plot[["samp_id"]])

    sashimi_plots <- list()

    for (i in seq_along(samp_ids)) {
        junctions_per_sample <- junctions_to_plot %>% dplyr::filter(samp_id == samp_ids[i])

        sashimi_plot <- gene_track_plot +
            ggplot2::geom_path(
                data = junctions_per_sample,
                ggplot2::aes(
                    x = x, y = y,
                    group = as.factor(index),
                    size = count,
                    colour = type
                ),
                lineend = "round",
                linetype = junctions_per_sample[["linetype"]]
            ) +
            ggplot2::scale_size_continuous(
                range = c(0.2, 1.5),
                limits = c(0, 1),
                guide = "none"
            ) +
            ggplot2::scale_linetype_manual(
                values = c(2, 1),
                guide = "none"
            ) +
            ggplot2::scale_colour_manual(
                name = "Junction type",
                values = annot_colour,
                breaks = c(
                    "annotated", "novel_acceptor", "novel_donor",
                    "novel_combo", "novel_exon_skip",
                    "ambig_gene", "unannotated"
                ),
                labels = c(
                    "Annotated", "Novel acceptor", "Novel donor",
                    "Novel combo", "Novel exon skip",
                    "Ambiguous gene", "Unannotated"
                ),
                drop = FALSE
            ) +
            ggplot2::ylab(samp_ids[i]) +
            ggplot2::theme(
                legend.key = ggplot2::element_rect(colour = "black", fill = "white"),
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_text(colour = "black", face = "bold", angle = 90)
            )

        if (count_label) {
            sashimi_plot <- sashimi_plot +
                ggrepel::geom_label_repel(
                    data = junctions_per_sample %>% dplyr::filter(mid_point),
                    ggplot2::aes(
                        x = x, y = y,
                        label = count,
                        colour = type
                    ),
                    min.segment.length = 0,
                    seed = 32,
                    show.legend = FALSE,
                    size = 3.5
                )
        }

        sashimi_plots[[i]] <- sashimi_plot
    }

    return(sashimi_plots)
}

#' Obtain coverage to plot
#'
#' `.coverage_to_plot_get` load and normalise the coverage from bigwigs to be
#' plotted and wrangle this into a format ready for `ggplot2`.
#'
#' @inheritParams plot_sashimi
#'
#' @param coords_to_plot list containing the coordinates to be used for
#'   plotting.
#'
#' @return A data.frame with the normalised coverage for cases and controls.
#'
#' @keywords internal
#' @noRd
.coverage_to_plot_get <- function(coords_to_plot,
    coverage_paths_case,
    coverage_paths_control,
    coverage_chr_control = NULL,
    load_func = .coverage_load,
    sum_func) {
    if (length(coverage_paths_case) > 1) {
        stop("Currently coverage plotting functions only support a single case sample as input")
    }

    pos <- coverage <- NULL

    # obtain coverage to plot on the x limits of plot
    region_to_plot <- GenomicRanges::GRanges(stringr::str_c(
        coords_to_plot[["chr"]], ":",
        coords_to_plot[["x_min"]], "-",
        coords_to_plot[["x_max"]]
    ))

    # obtain coverage for cases
    coverage_to_plot_case <- .coverage_load(
        coverage_path = coverage_paths_case,
        regions = region_to_plot,
        sum_fun = "mean",
        chr_format = NULL,
        method = "rt"
    ) %>%
        unlist() %>%
        dplyr::tibble(
            pos = start(region_to_plot):end(region_to_plot),
            coverage = .
        ) %>%
        dplyr::mutate(case_control = "case")

    # obtain coverage for controls
    if (!is.null(coverage_paths_control)) {
        coverage_to_plot_control <- vector(
            mode = "list",
            length = length(coverage_paths_control)
        )

        for (i in seq_along(coverage_paths_control)) {
            coverage_to_plot_control[[i]] <- .coverage_load(
                coverage_path = coverage_paths_control[i],
                regions = region_to_plot,
                sum_fun = "mean",
                chr_format = coverage_chr_control,
                method = "rt"
            ) %>%
                unlist() %>%
                dplyr::tibble(
                    pos = start(region_to_plot):end(region_to_plot),
                    coverage = .
                )
        }

        coverage_to_plot_control <- do.call(dplyr::bind_rows, coverage_to_plot_control) %>%
            dplyr::group_by(pos) %>%
            dplyr::summarise(coverage = sum_func(coverage)) %>%
            dplyr::mutate(case_control = "control")

        coverage_to_plot <- dplyr::bind_rows(
            coverage_to_plot_case,
            coverage_to_plot_control
        )
    } else {
        coverage_to_plot <- coverage_to_plot_case
    }

    # normalise coverage to a relative coverage across the region
    # to compare between case and controls
    coverage_to_plot <- coverage_to_plot %>%
        dplyr::group_by(case_control) %>%
        dplyr::mutate(coverage = coverage / sum(coverage)) %>%
        dplyr::ungroup()

    return(coverage_to_plot)
}

#' Plot coverage across region/transcript/gene of interest
#'
#' `.plot_coverage` will plot the coverage for case (and controls) across the
#' region of interest as a `ggplot`.
#'
#' @inheritParams plot_sashimi
#'
#' @param coverage_to_plot data.frame containing the normalised coverat to be
#'   plotted returned by `.coverage_to_plot_get`.
#' @param coords_to_plot list containing the coordinates to be used for
#'   plotting.
#'
#' @return An annotated sashimi plot.
#'
#' @keywords internal
#' @noRd
.plot_coverage <- function(coverage_to_plot, coords_to_plot, binwidth) {
    pos <- coverage <- NULL

    coverage_plot <-
        ggplot2::ggplot() +
        ggplot2::stat_summary_bin(
            data = coverage_to_plot,
            ggplot2::aes(
                x = pos,
                y = coverage,
                colour = case_control,
                fill = case_control
            ),
            fun = "mean",
            geom = "area",
            binwidth = binwidth,
            alpha = 0.2
        ) +
        ggplot2::scale_x_continuous(
            name = stringr::str_c("Chromosome ", coords_to_plot[["chr"]] %>%
                unique() %>%
                stringr::str_replace("chr", "")),
            limits = c(coords_to_plot[["x_min"]], coords_to_plot[["x_max"]])
        ) +
        ggplot2::scale_y_continuous(name = "Coverage") +
        ggplot2::scale_colour_manual(name = "", values = ggpubr::get_palette("jco", 2)) +
        ggplot2::scale_fill_manual(name = "", values = ggpubr::get_palette("jco", 2)) +
        ggpubr::theme_pubclean(flip = T) +
        ggplot2::theme(
            axis.line = ggplot2::element_line(colour = "black"),
            axis.ticks.x = ggplot2::element_blank()
        )

    return(coverage_plot)
}

#' Merge coverage and junction plots
#'
#' `.merge_coverage_sashimi` will combine the coverage and junction plots into
#' one. It will extract the legend from the sashimi plot in order to not
#' duplicate this for each junction plot.
#'
#' @param coverage_plot ggplot containing the coverage plot.
#' @param sashimi_plots list containing junction plots.
#'
#' @return sashimi plot.
#'
#' @keywords internal
#' @noRd
.merge_coverage_sashimi <- function(coverage_plot, sashimi_plots) {
    sashimi_legend <- ggpubr::get_legend(sashimi_plots[[1]]) %>%
        ggpubr::as_ggplot()

    sashimi_plots <- c(
        list(coverage_plot),
        sashimi_plots,
        list(sashimi_legend)
    )

    return(sashimi_plots)
}

#' Add annotation for sashimi plots
#'
#' `.plot_annotation` will add the gene name, chromosome and strand plotted.
#'
#' @inheritParams plot_sashimi
#'
#' @param sashimi_plots
#'   [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'    object containing junctions to be plotted.
#' @param coords_to_plot list containing the coordinates to be used for plotting.
#'
#' @return An annotated sashimi plot.
#'
#' @keywords internal
#' @noRd
.plot_annotation <- function(sashimi_plots, gene_tx_id, coords_to_plot) {
    sashimi_plots <- sashimi_plots %>%
        ggpubr::annotate_figure(
            top = ggpubr::text_grob(stringr::str_c(
                "Chromosome ",
                coords_to_plot[["chr"]] %>% stringr::str_replace("chr", ""),
                ", ", gene_tx_id,
                ", strand: ", coords_to_plot[["strand"]]
            ),
            x = 0.98, face = "italic", size = 10, just = "right"
            )
        )

    return(sashimi_plots)
}
