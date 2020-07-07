library(rdrop2)
library(tidyverse)
library(stringr)

##### Download and retrieve GTEx metadata from recount2 #####

gtex_metadata <-
    read_delim("http://snaptron.cs.jhu.edu/data/gtex/samples.tsv", delim = "\t")

gtex_metadata <- gtex_metadata %>%
    filter(
        SMTSD %in% c("Cells - Transformed fibroblasts"),
        SMAFRZE == "USE ME"
    ) %>%
    mutate(SMTSD_tidy = SMTSD %>%
        str_to_lower() %>%
        str_replace_all(" |-", "_") %>%
        str_replace_all("__*", "_"))

##### Download junction data using snaptron #####

# keep only chromosomes 1-22, X, Y and M
# fetch details/length of these
# chr_info <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC("hg38", ) %>%
#   dplyr::filter(UCSC_seqlevel %in% stringr::str_c(str_c("chr", c(1:22, "X", "Y", "M"))))

# server for the above not working on 2020/07/7
# so using pre downloaded chr_info
chr_info <- read_delim("/data/references/chr_info/hg38.chrom.sizes",
    delim = "\t",
    col_names = c("UCSC_seqlevel", "UCSC_seqlength")
) %>%
    dplyr::filter(UCSC_seqlevel %in% stringr::str_c(str_c("chr", c(1:22, "X", "Y", "M"))))

# downloads junction data for selected tissues using snaptron
tissues <- unique(gtex_metadata[["SMTSD_tidy"]])
GTEx_juncs_all <- list()

for (i in seq_along(tissues)) {
    GTEx_juncs <- tibble()

    gtex_metadata_tissue <- gtex_metadata %>%
        filter(SMTSD_tidy == tissues[i])

    # downloads 1 chromosome at a time
    for (j in 1:nrow(chr_info)) {
        print(str_c(Sys.time(), " - downloading junction data for ", chr_info[["UCSC_seqlevel"]][j]))

        juncs_per_chr <-
            read_delim(str_c(
                "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=",
                chr_info[["UCSC_seqlevel"]][j],
                ":1-",
                chr_info[["UCSC_seqlength"]][j],
                "&sids=", str_c(gtex_metadata_tissue[["rail_id"]], collapse = ",")
            ),
            delim = "\t"
            )

        GTEx_juncs <- bind_rows(GTEx_juncs, juncs_per_chr)
    }

    GTEx_juncs_all[[i]] <- GTEx_juncs
}

names(GTEx_juncs_all) <- tissues

##### Format and upload junction data #####

tidy_raw_count <- function(x) {
    x <- x[x != ""]
    y <- x %>%
        str_replace(".*:", "") %>%
        as.integer()

    names(y) <- x %>%
        str_replace(":.*", "")

    y <- y %>%
        t() %>%
        as_tibble()

    return(y)
}

# for now, upload these to dropbox to be able to loaded by dasper users
# potentially transfer to ExperimentalHub in future
drop_create("public/dasper/GTEx_v6_juncs")

for (i in 1:length(GTEx_juncs_all)) {

    # reformat GTEx junc counts into a count matrix
    GTEx_juncs_tidy <- GTEx_juncs_all[[i]][["samples"]] %>%
        str_split(pattern = ",")

    GTEx_juncs_tidy <- GTEx_juncs_tidy %>%
        lapply(FUN = tidy_raw_count) %>%
        do.call(args = ., bind_rows)

    # format colnames to be R-friendly
    colnames(GTEx_juncs_tidy) <- colnames(GTEx_juncs_tidy) %>% str_c("gtex_", .)

    # replace all NAs with 0s
    GTEx_juncs_tidy[is.na(GTEx_juncs_tidy)] <- 0

    # add back the original
    GTEx_juncs_tidy <- GTEx_juncs_tidy %>%
        bind_cols(GTEx_juncs_all[[i]] %>% dplyr::select(chr = chromosome, start, end, strand), .)

    # upload to dropbox
    save(GTEx_juncs_tidy,
        file = str_c("data/GTEx_juncs_", names(GTEx_juncs_all)[i], ".rda"),
        compress = "gzip"
    )

    drop_upload(str_c("data/GTEx_juncs_", names(GTEx_juncs_all)[i], ".rda"),
        path = "public/dasper/GTEx_v6_juncs"
    )

    file.remove(str_c("data/GTEx_juncs_", names(GTEx_juncs_all)[i], ".rda"))
}

# fibros - https://www.dropbox.com/s/6w3nrbzt3nknkmh/GTEx_juncs_cells_transformed_fibroblasts.rda?dl=0
