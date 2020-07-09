library(tidyverse)
library(stringr)

junctions <- read_delim("/data/RNA_seq_diag/mito/STAR/L1556-2624F_SJ.out.tab",
    delim = "\t",
    col_names = F,
    col_types = cols(X1 = "c")
)

# generate two sets of random 10,000 junctions for testing
for (i in 1:2) {
    set.seed(i)
    rand_indexes <- sort(sample(1:nrow(junctions), 10000))
    example_junctions <- junctions[rand_indexes, ]

    # store these in inst/extdata for example loading in the vignette
    # raw data should live here in according to http://r-pkgs.had.co.nz/data.html
    write_delim(example_junctions,
        str_c("inst/extdata/example_junctions_", i, ".txt"),
        delim = "\t",
        col_names = F
    )
}

# store .rda of example_junctions as in /data
# to be used by users
example_junctions_1_path <- system.file("extdata", "example_junctions_1.txt", package = "dasper", mustWork = TRUE)
example_junctions_2_path <- system.file("extdata", "example_junctions_2.txt", package = "dasper", mustWork = TRUE)

example_junctions <-
    junction_load(
        junction_paths = c(example_junctions_1_path, example_junctions_2_path),
        metadata = dplyr::tibble(samp_id = c("example_1", "example_2"))
    )

usethis::use_data(example_junctions, compress = "gzip")
