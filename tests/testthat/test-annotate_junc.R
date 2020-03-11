context("Annotating juncs with reference")

load("../../data/example_juncs.rda")
load("example_gtf.rda")

example_juncs_gr <-
  example_juncs %>%
  dplyr::select(chr, start, end, strand) %>%
  as.data.frame() %>%
  GRanges()

junc_metadata <- annotate_junc_ref(junc_metadata = example_juncs_gr, gtf = ref)

test_that("tidying annotation correctly infers strand", {

  x <- GRanges(c("1:1-100", "2:2-200", "3:3-300", "4:4-400"))
  strand(x) <- c("+", "-", "*", "*")
  mcols(x)["strand_start"] <- CharacterList("1" = "+", "2" = "*", "3" = "+", "4" = "*")
  mcols(x)["strand_end"] <- mcols(x)["strand_start"]

  expect_identical(as.character(strand(.tidy_junc_annot(x, cols_to_merge = "strand"))),
                   c("+", "-", "+", "*"))

})

test_that("the general junc_metadata output lookss correct", {
  expect_match(class(junc_metadata), "GRanges")
  expect_equal(length(junc_metadata), length(example_juncs_gr))
  expect_false(any(is.na(junc_metadata$junc_cat)))
})

test_that("annotate_junc_ref catches user-input errors", {
  expect_error(annotate_junc_ref(junc_metadata = example_juncs, gtf = ref),
               "junction_metadata must be in a GRanges format")
  expect_error(annotate_junc_ref(junc_metadata = example_juncs_gr, gtf = 40),
               "gtf must either be a path to the .gtf file or a pre-loaded gtf of class ensemblGenome")
})

test_that("exon annotation has been correctly retreived", {

  ref_exons <- ref@ev$gtf[ref@ev$gtf$feature == "exon",] %>%
    GRanges()

  # change random 20 junctions "manually" that transcript details match up to exon
  check <- TRUE

  for(i in sample(1:length(junc_metadata), 20)){

    junc_to_test <- junc_metadata[i]

    expect_exons_start <-
      ref_exons[end(ref_exons) == (start(junc_to_test) - 1) & as.character(seqnames(ref_exons)) == as.character(seqnames(junc_to_test))]
    expect_exons_end <-
      ref_exons[start(ref_exons) == (end(junc_to_test) + 1) & as.character(seqnames(ref_exons)) == as.character(seqnames(junc_to_test))]

    check <- all(check, identical(expect_exons_start$transcript_id %>% unique() %>% sort(),
                                  junc_to_test$transcript_id_start %>% unlist() %>% unname() %>% sort()))
    check <- all(check, identical(expect_exons_end$transcript_id %>% unique() %>% sort(),
                                  junc_to_test$transcript_id_end %>% unlist() %>% unname() %>% sort()))

  }

  expect_true(check)

})

test_that("junction categories meet expectations", {

  # split junctions by categories for
  annotated <- junc_metadata[junc_metadata$junc_cat == "annotated"]
  novel_donor <- junc_metadata[junc_metadata$junc_cat == "novel_donor"]
  novel_acceptor <- junc_metadata[junc_metadata$junc_cat == "novel_acceptor"]
  novel_exon_skip <- junc_metadata[junc_metadata$junc_cat == "novel_exon_skip"]
  novel_combo <- junc_metadata[junc_metadata$junc_cat == "novel_combo"]
  ambig_gene <- junc_metadata[junc_metadata$junc_cat == "ambig_gene"]
  none <- junc_metadata[junc_metadata$junc_cat == "none"]

  expect_true(all(lengths(annotated$gene_id_junc) > 0))
  expect_true(all(lengths(novel_donor[strand(novel_donor) == "+"]$gene_id_start) == 0))
  expect_true(all(lengths(novel_donor[strand(novel_donor) == "-"]$gene_id_end) == 0))
  expect_true(all(lengths(novel_acceptor[strand(novel_acceptor) == "+"]$gene_id_end) == 0))
  expect_true(all(lengths(novel_acceptor[strand(novel_acceptor) == "-"]$gene_id_start) == 0))
  expect_false(any(any(novel_combo$transcript_id_start %in% novel_combo$transcript_id_end)))
  expect_false(any(novel_exon_skip$junc_in_ref))
  expect_true(all(lengths(ambig_gene$gene_id_junc) > 1))
  expect_false(any(none$junc_in_ref))
  expect_true(length(c(unlist(none$exon_id_start), unlist(none$exon_id_end))) == 0)

})
