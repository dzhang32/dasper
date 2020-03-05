library(tidyverse)

juncs <- read_delim("/data/RNA_seq_diag/mito/STAR/L1556-2624F_SJ.out.tab",
                    delim = "\t",
                    col_names = F,
                    col_types = cols(X1 = "c"))

set.seed(32)
rand_indexes <- sort(sample(1:nrow(juncs), 10000))
example_juncs <- juncs[rand_indexes,]

write_delim(example_juncs, "data-raw/example_juncs.txt", delim = "\t", col_names = F)

colnames(example_juncs) <- c("chr", "start", "end", "strand", "intron_motif", "annotation", "uniq_map_read_count", "multi_map_read_count", "max_overhang")

usethis::use_data(example_juncs)
