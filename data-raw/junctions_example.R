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
    junctions_example <- junctions[rand_indexes, ]

    # store these in inst/extdata for example loading in the vignette
    # raw data should live here in according to http://r-pkgs.had.co.nz/data.html
    write_delim(junctions_example,
        str_c("inst/extdata/junctions_example_", i, ".txt"),
        delim = "\t",
        col_names = F
    )
}

# generate a junctions_example.rda for testing/vignette and as illustration of expected junction format for users
junctions_example_1_path <- system.file("extdata", "junctions_example_1.txt", package = "dasper", mustWork = TRUE)
junctions_example_2_path <- system.file("extdata", "junctions_example_2.txt", package = "dasper", mustWork = TRUE)

# only take chr21 + 22
junctions_example <-
    junction_load(
        junction_paths = c(junctions_example_1_path, junctions_example_2_path),
        controls = "fibroblasts",
        chrs = c("21", "22")
    )

# and first 3 GTEx control samples to save space
raw_counts <- junctions_example[, c(1:5)]

usethis::use_data(junctions_example, compress = "xz")
