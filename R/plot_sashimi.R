#' Visualise RNA-seq data in a the form of a sashimi plot
#'
#' \code{plot_sashimi} plots the splicing events and coverage over a specific
#' genes/transcript and/or genomic region of interest. The plotting are built
#' from `ggplot2` functions.
#'
#' @inheritParams junction_annot
#'
#' @param gene_tx_id gene name as ensembl id, ensembl transcript id or gene
#'   name.
#' @param case_id list of one element. This must be a character vector
#'   containing the
#' @param control_agg_func list of one element. This must be a character vector
#'   containing the
#' @param region a [GenomicRanges][GenomicRanges::GRanges-class] of length 1
#'   that is used to filter the exons/junctions. Only those that overlap this
#'   region are plotted.
#' @param annot_colour character vector of colours for junction types. One value
#'   must be supplied labelling each of the 7 possible types.
#' @param digits used in `round`, specifying the number of digits to round the
#'   junction counts to for visualisation purposes.
#' @param count_label logical value specifying whether to add label the count of
#'   each junction.
#'
#' @return `ggplot` displaying the splicing (and coverage) surrounding the
#'   transcript/region of interest.
#'
#' @family plotting
#' @export
#'
#' @examples
#'
#' if (!exists("ref")) {
#'     # use Genomic state to load txdb (GENCODE v31)
#'     ref <- GenomicState::GenomicStateHub(version = "31", genome = "hg38", filetype = "TxDb")[[1]]
#'     # convert seqlevels to match junctions
#'     seqlevels(ref) <- stringr::str_replace(seqlevels(ref), "chr", "")
#' }
#'
#' if (!exists("junctions_processed")) {
#'     junctions_processed <-
#'         junction_process(
#'             junctions_example,
#'             ref,
#'             count_thresh = c("raw" = 5),
#'             n_samp = c("raw" = 1),
#'             width_range = c(25, 1000000),
#'             types = c("ambig_gene", "unannotated"),
#'         )
#'     GenomeInfoDb::seqlevels(junctions_processed) <-
#'         paste0("chr", GenomeInfoDb::seqlevels(junctions_processed))
#' }
#'
#' sashimi_plot <- plot_sashimi(
#'     junctions = junctions_processed,
#'     ref = ref,
#'     gene_tx_id = "ENSG00000241973.10"
#' )
plot_sashimi <- function(junctions,
    ref,
    gene_tx_id,
    case_id = NULL,
    control_agg_func = mean,
    region = NULL,
    annot_colour = c(
        ggpubr::get_palette("jco", 1),
        ggpubr::get_palette("npg", 7)[c(1, 3, 2, 5, 6)],
        ggpubr::get_palette("jco", 6)[c(3)]
    ),
    digits = 2,
    count_label = TRUE) {

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

    sashimi_plots <- .plot_sashimi_junctions(junctions_to_plot,
        gene_track_plot,
        case_id,
        control_agg_func,
        digits,
        assay_name = "norm",
        annot_colour,
        count_label
    )

    ##### Plot coverage #####

    # if(!is.null(coverage_case)){
    #
    #   coverage_plot <-
    #     plot_coverage(coverage_case = coverage_case,
    #                   coverage_ctrl = coverage_ctrl,
    #                   region_to_plot = GenomicRanges::GRanges(stringr::str_c(chromosome_to_plot, ":", xmin, "-", xmax)),
    #                   binwidth = binwidth) +
    #     theme(axis.title.x = element_blank())
    #
    #   coverage_plot_leg <- ggpubr::get_legend(coverage_plot)
    #   sashimi_plot_plot_leg <- ggpubr::get_legend(sashimi_plot_list[[1]])
    #
    #   # arrange plots into panels
    #   sashimi_plot <- ggpubr::ggarrange(plotlist = c(list(coverage_plot), sashimi_plot_list),
    #                                     ncol = 1,
    #                                     nrow = length(c(list(coverage_plot), sashimi_plot_list)),
    #                                     align = "v",
    #                                     legend = "none",
    #                                     heights = c(1, rep(1.25, length(sashimi_plot_list))))
    #
    #   # add legends
    #   sashimi_plot <- ggpubr::ggarrange(plotlist = list(coverage_plot_leg, sashimi_plot, sashimi_plot_plot_leg),
    #                                     ncol = 1,
    #                                     nrow = 3,
    #                                     heights = c(0.75, 10, 0.75))
    #
    # }

    ##### Add annotation #####

    sashimi_plots <- .plot_annotation(sashimi_plots, gene_tx_id, coords_to_plot)

    return(sashimi_plots)
}

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

.gene_tx_to_plot_get <- function(ref,
    gene_tx_id,
    gene_tx_type) {

    # create named list of gene/tx id
    gene_tx_id_list <- list(gene_tx_id)
    names(gene_tx_id_list) <- gene_tx_type

    # filter for exons of gene/tx of interest
    gene_tx_to_plot <- GenomicFeatures::genes(ref, filter = gene_tx_id_list)

    return(gene_tx_to_plot)
}

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

.junctions_to_plot_get <- function(junctions, gene_tx_list, region) {
    gene_tx <- gene_tx_list %>% unlist()

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
    coords_to_plot[["x_min"]] <- coords_to_plot[["min_exon_start"]] - coords_to_plot[["range_exon_start_end"]] / ext_factor
    coords_to_plot[["x_max"]] <- coords_to_plot[["max_exon_end"]] + coords_to_plot[["range_exon_start_end"]] / ext_factor

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

    return(coords_to_plot)
}


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
            limits = c(
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

.plot_sashimi_junctions <- function(junctions_to_plot,
    gene_track_plot,
    case_id,
    control_agg_func,
    digits,
    assay_name,
    annot_colour,
    count_label) {

    # format junctions into df with count/type details
    junctions_to_plot <- .junctions_counts_type_get(
        junctions_to_plot = junctions_to_plot,
        case_id = case_id,
        control_agg_func = control_agg_func,
        digits = digits,
        assay_name = assay_name
    )

    # obtain points for curves of junctions
    junctions_to_plot <- .junctions_points_get(junctions_to_plot, ncp = 25)

    sashimi_plots <- .plot_junctions(
        junctions_to_plot = junctions_to_plot,
        gene_track_plot = gene_track_plot,
        annot_colour = annot_colour,
        count_label = count_label
    )

    sashimi_plots <- ggpubr::ggarrange(
        plotlist = sashimi_plots,
        ncol = 1,
        nrow = length(sashimi_plots),
        common.legend = TRUE,
        align = "v",
        legend = "bottom"
    )
}

.junctions_counts_type_get <- function(junctions_to_plot, case_id = list(samp_id = "samp_1"), control_agg_func = mean, digits = 2, assay_name = "norm") {

    # for R CMD Check
    index <- type <- . <- NULL

    junctions_counts <- junctions_to_plot %>%
        GenomicRanges::ranges() %>%
        as.data.frame() %>%
        dplyr::mutate(
            index = dplyr::row_number(),
            type = mcols(junctions_to_plot)[["type"]]
        )

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
    if (!is.null(control_agg_func)) {
        which_control <- which(colData(junctions_to_plot)[["case_control"]] == "control")

        junctions_counts[["control"]] <-
            assays(junctions_to_plot)[[assay_name]][, which_control] %>%
            apply(MARGIN = 1, FUN = control_agg_func) %>%
            round(digits = digits)
    }

    # filter for out junctions that are not expressed (> 0 counts) in any sample
    junctions_counts <- junctions_counts %>%
        dplyr::select(-start, -end, -width, -index, -type) %>%
        apply(MARGIN = 1, FUN = function(x) !all(x == 0)) %>%
        dplyr::filter(junctions_counts, .)

    return(junctions_counts)
}

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

.plot_annotation <- function(sashimi_plots, gene_tx_id, coords_to_plot) {
    sashimi_plots <- sashimi_plots %>%
        ggpubr::annotate_figure(
            top = ggpubr::text_grob(stringr::str_c(
                "Chromosome ", coords_to_plot[["chr"]], ", ",
                gene_tx_id, ", strand: ", coords_to_plot[["strand"]]
            ),
            x = 0.98, face = "italic", size = 10, just = "right"
            )
        )

    return(sashimi_plots)
}
