library(tidyverse)

# for testing, we require a gtf loaded in through refGenome for annotate_junc_ref
gtf <- "/data/references/ensembl/gtf_gff3/v95/Homo_sapiens.GRCh38.95.gtf"
ref <- refGenome::ensemblGenome()
refGenome::basedir(ref) <- dirname(gtf)
refGenome::read.gtf(ref, gtf %>% stringr::str_replace(".*/", ""))

# place this in testthat as does not need to be exported with package
save(ref, file = "tests/testthat/example_gtf.rda", compress = "gzip")

# check compressed size/type
# tools::checkRdaFiles("tests/testthat/example_gtf.rda")
