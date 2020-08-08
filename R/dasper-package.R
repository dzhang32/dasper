#' dasper: detecting abberant splicing events from RNA-seq data
#'
#' Placeholder for package description - to be updated
#'
#' @docType package
#' @name dasper
NULL

#' @import S4Vectors
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @importFrom dplyr count
#' @importFrom GenomicRanges seqnames start end strand width
#' @importFrom GenomicRanges start<- end<- strand<-
#' @importFrom GenomicRanges findOverlaps GRanges
#' @importFrom IRanges CharacterList
#' @importFrom SummarizedExperiment colData assays
#' @importFrom SummarizedExperiment assays<-
NULL

##### Using sklearn in R via reticulate #####

# global reference to scipy (will be initialized in .onLoad)
sklearn <- NULL

#' Load sklearn upon initialisation of dasper
#'
#' \code{.onLoad} uses superassignment to update global reference to sklearn.
#' Delay load means allows successfull load without sklearn installed. This
#' allows users to set python environment using reticulate::use_virtualenv().
#'
#' @param libname placeholder
#' @param pkgname placeholder
#'
#' @seealso https://rstudio.github.io/reticulate/articles/package.html
.onLoad <- function(libname, pkgname) {
    sklearn <<- reticulate::import("sklearn",
        delay_load = TRUE
    )
}
