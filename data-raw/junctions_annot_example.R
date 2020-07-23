library(tidyverse)
library(stringr)

junctions_annot_example <- junction_annot(junctions_example,
    ref = "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
)

usethis::use_data(junctions_annot_example, compress = "xz")
