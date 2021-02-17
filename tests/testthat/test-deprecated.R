# The below code is used to test the deprecated function annot_junc_ref()
# to save time, this is not run as part of the automatic testing of the package
#
# junctions_gr <- junctions_example %>% rowRanges()
# GenomeInfoDb::seqlevels(junctions_gr) <- GenomeInfoDb::seqlevels(junctions_gr) %>% stringr::str_remove("chr")
#
# junctions_annot_old <- annotate_junc_ref(junctions_gr, "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf")
# junctions_annot_new <- junction_annot(junctions_gr, "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf")
#
# junctions_annot_old$junc_cat %>% table()
# junctions_annot_new$type %>% table()
