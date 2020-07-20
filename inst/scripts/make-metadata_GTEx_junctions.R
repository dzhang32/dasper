library(sessioninfo)
library(stringr)

tissues <- c("fibroblast", "whole_blood", "skeletal_muscle", "lymphocytes", "skin_sun_exposed")

gtex_junctions_metadata <- data.frame(
  Title = str_c("GTEx_v6_junctions_",
                tissues),
  Description = str_c(
    "data.frame containing counts for all junctions originating from GTEx v6 ", str_replace_all(tissues, "_", " "), " samples.",
    "Each row corresponds to a junction. Columns detail the 'chr', 'start', 'end', 'strand' and remaining columns refer to each GTEx sample in the format 'gtex_XXXXX', XXXXX being the rail_id from snaptron.",
    "Data has been released via the recount2 project, downloaded by snaptron and wrangled into convienient format for use in the dasper package."
  ),
  BiocVersion = "3.11",
  Genome = "GRCh38",
  SourceType = "TSV",
  SourceUrl = "http://snaptron.cs.jhu.edu/",
  SourceVersion = "Jul 20 2020",
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = "Genotype-Tissue Expression (GTEx) project - https://www.gtexportal.org/home/",
  Maintainer = "David Zhang <dyzhang32@gmail.com>",
  RDataClass = c("data.frame"),
  DispatchClass = "Rda",
  RDataPath = file.path("dasper",
                        str_c("GTEx_v6_junctions_",
                              tissues,
                              ".rda")
                        ),
  Tags = "GTEx_v6_junctions",
  row.names = NULL,
  stringsAsFactors = FALSE
)

write.csv(gtex_junctions_metadata, file = "inst/extdata/metadata.csv", row.names = FALSE)

# check metadata format complies with Bioconductor guidelines
# AnnotationHubData::makeAnnotationHubMetadata(getwd())

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 3.6.1 (2019-07-05)
# os       Ubuntu 16.04.6 LTS
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_GB.UTF-8
# ctype    en_GB.UTF-8
# tz       Europe/London
# date     2020-07-20
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version   date       lib source
# AnnotationDbi            1.46.1    2019-08-20 [2] Bioconductor
# AnnotationForge          1.26.0    2019-05-02 [1] Bioconductor
# AnnotationHub            2.16.1    2019-09-04 [1] Bioconductor
# AnnotationHubData        1.14.0    2019-05-02 [1] Bioconductor
# assertthat               0.2.1     2019-03-21 [2] CRAN (R 3.6.1)
# Biobase                  2.44.0    2019-05-02 [2] Bioconductor
# BiocFileCache            1.8.0     2019-05-02 [1] Bioconductor
# BiocGenerics             0.30.0    2019-05-02 [2] Bioconductor
# BiocManager              1.30.10   2019-11-16 [1] CRAN (R 3.6.1)
# BiocParallel             1.18.1    2019-08-06 [1] Bioconductor
# biocViews                1.52.2    2019-05-29 [1] Bioconductor
# biomaRt                  2.40.5    2019-10-01 [1] Bioconductor
# Biostrings               2.52.0    2019-05-02 [2] Bioconductor
# bit                      1.1-14    2018-05-29 [2] CRAN (R 3.6.1)
# bit64                    0.9-7     2017-05-08 [2] CRAN (R 3.6.1)
# bitops                   1.0-6     2013-08-17 [2] CRAN (R 3.6.1)
# blob                     1.2.0     2019-07-09 [2] CRAN (R 3.6.1)
# cli                      2.0.2     2020-02-28 [1] CRAN (R 3.6.1)
# crayon                   1.3.4     2017-09-16 [2] CRAN (R 3.6.1)
# curl                     4.3       2019-12-02 [1] CRAN (R 3.6.1)
# data.table               1.12.8    2019-12-09 [1] CRAN (R 3.6.1)
# DBI                      1.1.0     2019-12-15 [2] CRAN (R 3.6.1)
# dbplyr                   1.4.2     2019-06-17 [2] CRAN (R 3.6.1)
# DelayedArray             0.10.0    2019-05-02 [2] Bioconductor
# digest                   0.6.25    2020-02-23 [1] CRAN (R 3.6.1)
# dplyr                    1.0.0     2020-05-29 [1] CRAN (R 3.6.1)
# ellipsis                 0.3.1     2020-05-15 [1] CRAN (R 3.6.1)
# fansi                    0.4.1     2020-01-08 [1] CRAN (R 3.6.1)
# fastmap                  1.0.1     2019-10-08 [1] CRAN (R 3.6.1)
# formatR                  1.7       2019-06-11 [2] CRAN (R 3.6.1)
# futile.logger          * 1.4.3     2016-07-10 [2] CRAN (R 3.6.1)
# futile.options           1.0.1     2018-04-20 [2] CRAN (R 3.6.1)
# generics                 0.0.2     2018-11-29 [1] CRAN (R 3.6.1)
# GenomeInfoDb             1.20.0    2019-05-02 [1] Bioconductor
# GenomeInfoDbData         1.2.1     2019-09-21 [2] Bioconductor
# GenomicAlignments        1.20.1    2019-06-18 [2] Bioconductor
# GenomicFeatures          1.36.4    2019-07-10 [1] Bioconductor
# GenomicRanges            1.36.1    2019-09-06 [2] Bioconductor
# glue                     1.4.1     2020-05-13 [1] CRAN (R 3.6.1)
# graph                    1.62.0    2019-05-02 [1] Bioconductor
# hms                      0.5.3     2020-01-08 [1] CRAN (R 3.6.1)
# htmltools                0.5.0     2020-06-16 [1] CRAN (R 3.6.1)
# httpuv                   1.5.2     2019-09-11 [2] CRAN (R 3.6.1)
# httr                     1.4.1     2019-08-05 [2] CRAN (R 3.6.1)
# interactiveDisplayBase   1.22.0    2019-05-02 [1] Bioconductor
# IRanges                  2.18.3    2019-09-24 [1] Bioconductor
# jsonlite                 1.7.0     2020-06-25 [1] CRAN (R 3.6.1)
# lambda.r                 1.2.4     2019-09-18 [2] CRAN (R 3.6.1)
# later                    1.1.0.1   2020-06-05 [1] CRAN (R 3.6.1)
# lattice                  0.20-38   2018-11-04 [4] CRAN (R 3.6.1)
# lifecycle                0.2.0     2020-03-06 [1] CRAN (R 3.6.1)
# magrittr                 1.5       2014-11-22 [2] CRAN (R 3.6.1)
# Matrix                   1.2-18    2019-11-27 [1] CRAN (R 3.6.1)
# matrixStats              0.55.0    2019-09-07 [2] CRAN (R 3.6.1)
# memoise                  1.1.0     2017-04-21 [2] CRAN (R 3.6.1)
# mime                     0.9       2020-02-04 [1] CRAN (R 3.6.1)
# OrganismDbi              1.26.0    2019-05-02 [1] Bioconductor
# pillar                   1.4.4     2020-05-05 [1] CRAN (R 3.6.1)
# pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 3.6.1)
# prettyunits              1.1.1     2020-01-24 [1] CRAN (R 3.6.1)
# progress                 1.2.2     2019-05-16 [2] CRAN (R 3.6.1)
# promises                 1.1.1     2020-06-09 [1] CRAN (R 3.6.1)
# purrr                    0.3.4     2020-04-17 [1] CRAN (R 3.6.1)
# R6                       2.4.1     2019-11-12 [1] CRAN (R 3.6.1)
# rappdirs                 0.3.1     2016-03-28 [1] CRAN (R 3.6.1)
# RBGL                     1.60.0    2019-05-02 [1] Bioconductor
# rBiopaxParser            2.24.0    2019-05-02 [1] Bioconductor
# Rcpp                     1.0.4.6   2020-04-09 [1] CRAN (R 3.6.1)
# RCurl                    1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)
# rlang                    0.4.6     2020-05-02 [1] CRAN (R 3.6.1)
# Rsamtools                2.0.3     2019-10-10 [1] Bioconductor
# RSQLite                  2.1.2     2019-07-24 [2] CRAN (R 3.6.1)
# rstudioapi               0.11      2020-02-07 [2] CRAN (R 3.6.1)
# rtracklayer              1.44.4    2019-09-06 [2] Bioconductor
# RUnit                    0.4.32    2018-05-18 [1] CRAN (R 3.6.1)
# S4Vectors                0.22.1    2019-09-09 [2] Bioconductor
# sessioninfo            * 1.1.1     2018-11-05 [2] CRAN (R 3.6.1)
# shiny                    1.5.0     2020-06-23 [1] CRAN (R 3.6.1)
# stringi                  1.4.6     2020-02-17 [1] CRAN (R 3.6.1)
# stringr                  1.4.0     2019-02-10 [2] CRAN (R 3.6.1)
# SummarizedExperiment     1.14.1    2019-07-31 [1] Bioconductor
# tibble                   3.0.1     2020-04-20 [1] CRAN (R 3.6.1)
# tidyselect               1.1.0     2020-05-11 [1] CRAN (R 3.6.1)
# vctrs                    0.3.1     2020-06-05 [1] CRAN (R 3.6.1)
# withr                    2.2.0     2020-04-20 [1] CRAN (R 3.6.1)
# XML                      3.99-0.3  2020-01-20 [2] CRAN (R 3.6.1)
# xtable                   1.8-4     2019-04-21 [1] CRAN (R 3.6.1)
# XVector                  0.24.0    2019-05-02 [2] Bioconductor
# yaml                     2.2.1     2020-02-01 [1] CRAN (R 3.6.1)
# zlibbioc                 1.30.0    2019-05-02 [2] Bioconductor
#
# [1] /home/dzhang/R/x86_64-pc-linux-gnu-library/3.6
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
