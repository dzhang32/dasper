library(tidyverse)

##### GTF - testing annotate_junc_ref #####

# for testing, we require a gtf loaded in through refGenome for annotate_junc_ref
gtf <- "/data/references/ensembl/gtf_gff3/v95/Homo_sapiens.GRCh38.95.gtf"
ref <- refGenome::ensemblGenome()
refGenome::basedir(ref) <- dirname(gtf)
refGenome::read.gtf(ref, gtf %>% stringr::str_replace(".*/", ""))

# place this in testthat as does not need to be exported with package
save(ref, file = "tests/testthat/example_gtf.rda", compress = "gzip")

##### Annotated juncs - testing filter_junc #####

load("data/example_juncs.rda")

example_juncs[["metadata"]] <-
  annotate_junc_ref(junc_metadata = example_juncs[["metadata"]],
                    gtf = ref)

# add "w_annot" to suffix to differentiate from example_juncs
example_juncs_w_annot <- example_juncs

save(example_juncs_w_annot, file = "tests/testthat/example_juncs_w_annot.rda", compress = "gzip")
