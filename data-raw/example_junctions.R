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
