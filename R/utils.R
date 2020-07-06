#' Filter a df for chromosomes to keep
#'
#' \code{.chr_filter} will filter a df, keeping only chromosomes listed in the
#' chrs argument. It includes a couple of checks to notify users whether some or
#' all chrs are not found in their data. This aims to catch situations where
#' alternative nomenclature may be used - e.g. "M" vs "MT", or "chr1" vs "1".
.chr_filter <- function(x, chrs) {

    # check if different chromosome formats
    if ((!any(chrs %in% x[["chr"]]))) {
        stop("No chromosomes in chrs match those of the junction data.")
    } else if (!(all(chrs %in% x[["chr"]]))) {
        warning(stringr::str_c(
            "The following chromosomes in chrs are not found in the junction data: ",
            stringr::str_c(chrs[!(chrs %in% x[["chr"]])], collapse = ", ")
        ))
    }

    x <- x %>% dplyr::filter(chr %in% chrs)

    return(x)
}
