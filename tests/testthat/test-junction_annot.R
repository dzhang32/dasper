context("Testing the annotation of junctions")

##### .junction_annot_tidy #####

x <- GRanges(c("1:1-100", "2:2-200", "3:3-300", "4:4-400"))
strand(x) <- c("+", "-", "*", "*")
mcols(x)["strand_start"] <- CharacterList("1" = "+", "2" = "*", "3" = "+", "4" = "*")
mcols(x)["strand_end"] <- mcols(x)["strand_start"]
strand_tidy <- as.character(strand(.junction_annot_tidy(x, cols_to_merge = "strand")))

test_that(".junction_annot_tidy correctly infers strand", {
    expect_identical(strand_tidy, c("+", "-", "+", "*"))
})

##### junction_annot #####

suppressWarnings(expr = {
    ref <- GenomicFeatures::makeTxDbFromGFF("ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz")
})

junctions <- junction_annot(junctions_example, ref)


test_that("junction_annot output generally looks correct", {
    expect_true(methods::isClass(junctions, "RangedSummarisedExperiment"))
    expect_false(any(is.na(mcols(junctions)[["type"]])))
})

test_that("junction_annot catches user-input errors", {
    expect_error(
        junction_annot(junctions = "not_a_RSE", ref = ref),
        "junctions must be in a RangedSummarisedExperiment format"
    )
    expect_error(
        junction_annot(junctions = junctions, ref = 50),
        "ref must either be a path to the .gtf/gff3 file or a pre-loaded TxDb object"
    )
})

# split junctions by categories
annotated <- junctions[mcols(junctions)[["type"]] == "annotated"]
novel_donor <- junctions[mcols(junctions)[["type"]] == "novel_donor"]
novel_acceptor <- junctions[mcols(junctions)[["type"]] == "novel_acceptor"]
novel_exon_skip <- junctions[mcols(junctions)[["type"]] == "novel_exon_skip"]
novel_combo <- junctions[mcols(junctions)[["type"]] == "novel_combo"]
ambig_gene <- junctions[mcols(junctions)[["type"]] == "ambig_gene"]
unannotated <- junctions[mcols(junctions)[["type"]] == "unannotated"]

test_that("junction categories meet expectations", {
    expect_identical(
        levels(mcols(junctions)[["type"]]),
        c(
            "annotated",
            "novel_acceptor",
            "novel_donor",
            "novel_exon_skip",
            "novel_combo",
            "ambig_gene",
            "unannotated"
        )
    )

    # annotated
    expect_true(all(lengths(mcols(annotated)[["gene_id_junction"]]) > 0))
    expect_true(all(lengths(mcols(annotated)[["gene_id_start"]]) > 0 &
        lengths(mcols(annotated)[["gene_id_end"]]) > 0))

    # novel_donor
    expect_true(all(lengths(mcols(novel_donor[strand(novel_donor) == "+"])[["gene_id_start"]]) == 0))
    expect_true(all(lengths(mcols(novel_donor[strand(novel_donor) == "+"])[["gene_id_end"]]) > 0))
    expect_true(all(lengths(mcols(novel_donor[strand(novel_donor) == "-"])[["gene_id_start"]]) > 0))
    expect_true(all(lengths(mcols(novel_donor[strand(novel_donor) == "-"])[["gene_id_end"]]) == 0))

    # novel_acceptor
    expect_true(all(lengths(mcols(novel_acceptor[strand(novel_acceptor) == "+"])[["gene_id_start"]]) > 0))
    expect_true(all(lengths(mcols(novel_acceptor[strand(novel_acceptor) == "+"])[["gene_id_end"]]) == 0))
    expect_true(all(lengths(mcols(novel_acceptor[strand(novel_acceptor) == "-"])[["gene_id_start"]]) == 0))
    expect_true(all(lengths(mcols(novel_acceptor[strand(novel_acceptor) == "-"])[["gene_id_end"]]) > 0))

    # novel_combo
    expect_false(any(any(mcols(novel_combo)[["tx_name_start"]] %in% mcols(novel_combo)[["tx_name_end"]])))

    # novel_exon_skip
    expect_true(all(any(mcols(novel_exon_skip)[["tx_name_start"]] %in% mcols(novel_exon_skip)[["tx_name_end"]])))

    # ambig_gene
    expect_true(all(lengths(ambig_gene$gene_id_junc) > 1))

    # unannotated
    expect_false(any(mcols(unannotated)[["in_ref"]]))
    expect_true(length(c(
        unlist(mcols(unannotated)[["exon_name_start"]]),
        unlist(mcols(unannotated)[["exon_name_end"]])
    )) == 0)
})

ref_exons <- ref %>% GenomicFeatures::exons(columns = c("gene_id", "tx_name", "exon_name"))
ref_exons_start_end <- .get_start_end(ref_exons)

# define function to manually check n randomly sampled junctions
# to make exon annotation from junction_annot matches expections
annot_check <- function(junctions, ref_exons, n) {
    check <- TRUE

    start(junctions) <- start(junctions) - 1
    end(junctions) <- end(junctions) + 1

    junctions_start_end <- .get_start_end(junctions)
    ref_exons_start_end <- .get_start_end(ref_exons)

    for (i in sample(1:length(junctions), n, replace = T)) {
        junction_to_test <- junctions_start_end %>%
            lapply(FUN = function(x) {
                x[i]
            })

        expect_exons_start <-
            ref_exons[findOverlaps(ref_exons_start_end[["end"]],
                junction_to_test[["start"]],
                ignore.strand = F
            ) %>%
                queryHits()]

        expect_exons_end <-
            ref_exons[findOverlaps(ref_exons_start_end[["start"]],
                junction_to_test[["end"]],
                ignore.strand = F
            ) %>%
                queryHits()]

        # check exons
        check <- all(check, identical(
            expect_exons_start$exon_name %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["exon_name_start"]] %>% unlist() %>% unname() %>% sort()
        ))

        check <- all(check, identical(
            expect_exons_end$exon_name %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["exon_name_end"]] %>% unlist() %>% unname() %>% sort()
        ))

        # check transcripts
        check <- all(check, identical(
            expect_exons_start$tx_name %>% unlist() %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["tx_name_start"]] %>% unlist() %>% unname() %>% sort()
        ))

        check <- all(check, identical(
            expect_exons_end$tx_name %>% unlist() %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["tx_name_end"]] %>% unlist() %>% unname() %>% sort()
        ))

        # check genes
        check <- all(check, identical(
            expect_exons_start$gene_id %>% unlist() %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["gene_id_start"]] %>% unlist() %>% unname() %>% sort()
        ))

        check <- all(check, identical(
            expect_exons_end$gene_id %>% unlist() %>% unique() %>% sort(),
            mcols(junction_to_test[["start"]])[["gene_id_end"]] %>% unlist() %>% unname() %>% sort()
        ))
    }

    return(check)
}

test_that("exon annotation has been correctly retreived", {
    expect_true(annot_check(annotated, ref_exons, 10))
    expect_true(annot_check(novel_donor, ref_exons, 10))
    expect_true(annot_check(novel_acceptor, ref_exons, 10))
    expect_true(annot_check(novel_exon_skip, ref_exons, 10))
    expect_true(annot_check(novel_combo, ref_exons, 5))
    expect_true(annot_check(ambig_gene, ref_exons, 5))
    expect_true(annot_check(unannotated, ref_exons, 10))
})
