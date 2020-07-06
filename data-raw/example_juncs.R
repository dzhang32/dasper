library(tidyverse)
library(stringr)

juncs <- read_delim("/data/RNA_seq_diag/mito/STAR/L1556-2624F_SJ.out.tab",
    delim = "\t",
    col_names = F,
    col_types = cols(X1 = "c")
)

# generate two sets of random 10,000 juncs for testing
for (i in 1:2) {
    set.seed(i)
    rand_indexes <- sort(sample(1:nrow(juncs), 10000))
    example_juncs <- juncs[rand_indexes, ]

    # store these in inst/extdata for example loading in the vignette
    # raw data should live here in according to http://r-pkgs.had.co.nz/data.html
    write_delim(example_juncs,
        str_c("inst/extdata/example_juncs_", i, ".txt"),
        delim = "\t",
        col_names = F
    )
}

example_juncs_1_path <- system.file("extdata", "example_juncs_1.txt", package = "dasper", mustWork = TRUE)
example_juncs_2_path <- system.file("extdata", "example_juncs_2.txt", package = "dasper", mustWork = TRUE)

example_juncs <-
    junc_load(
        junc_paths = c(example_juncs_1_path, example_juncs_2_path),
        metadata = dplyr::tibble(samp_id = c("example_1", "example_2"))
    )

usethis::use_data(example_juncs, compress = "gzip")
