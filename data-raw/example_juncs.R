library(tidyverse)
library(stringr)

juncs <- read_delim("/data/RNA_seq_diag/mito/STAR/L1556-2624F_SJ.out.tab",
                    delim = "\t",
                    col_names = F,
                    col_types = cols(X1 = "c"))

# generate two sets of random 10,000 juncs for testing
for(i in 1:2){

  set.seed(i)
  rand_indexes <- sort(sample(1:nrow(juncs), 10000))
  example_juncs <- juncs[rand_indexes,]

  write_delim(example_juncs,
              str_c("data-raw/example_juncs_", i, ".txt"),
              delim = "\t",
              col_names = F)

}

# merged juncs into /data for vignette
example_juncs <-
  dasper::merge_junc(junc_paths = c("data-raw/example_juncs_1.txt",
                                    "data-raw/example_juncs_2.txt"),
                     sample_ids = c("eg1", "eg2"))

usethis::use_data(example_juncs)
