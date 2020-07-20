library(rdrop2)
library(tidyverse)
library(stringr)

##### Download and retrieve GTEx metadata from recount2 #####

# generate df detailing key between gtex name and tidy name
# for all tissues of interest (CATs)
tissues_of_interest <-
  tibble(tissue_tidy = c("fibroblast",
                         "whole_blood",
                         "skeletal_muscle",
                         "lymphocytes",
                         "skin_sun_exposed"),
         tissue_gtex = c("Cells - Transformed fibroblasts",
                         "Whole Blood",
                         "Muscle - Skeletal",
                         "Cells - EBV-transformed lymphocytes",
                         "Skin - Sun Exposed (Lower leg)"))

# download and filter GTEx metadata
gtex_metadata <- read_delim("http://snaptron.cs.jhu.edu/data/gtex/samples.tsv", delim = "\t")

if(!all(tissues_of_interest[["tissue_gtex"]] %in% gtex_metadata[["SMTSD"]])){

  stop("Not all tissues of interest found in GTEx metadata")

}

gtex_metadata <- gtex_metadata %>%
    filter(
        SMTSD %in% tissues_of_interest[["tissue_gtex"]],
        SMAFRZE == "USE ME"
    )

##### Download and tidy junction data using snaptron #####

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
# drop_auth()
drop_create("public/dasper/GTEx_v6_junctions")

# downloads junction data for selected tissues using snaptron
for (i in 1:nrow(tissues_of_interest)) {
    GTEx_junctions <- tibble()

    print(str_c(Sys.time(), " - ", tissues_of_interest[["tissue_tidy"]][i]))

    gtex_metadata_tissue <- gtex_metadata %>%
        filter(SMTSD == tissues_of_interest[["tissue_gtex"]][i])

    # downloads junctions 1 chromosome at a time
    for (j in 1:nrow(chr_info)) {
        print(str_c(Sys.time(), " - downloading junction data for ", chr_info[["UCSC_seqlevel"]][j]))

      suppressMessages(expr = {
        junctions_per_chr <-
          read_delim(str_c(
            "http://snaptron.cs.jhu.edu/gtex/snaptron?regions=",
            chr_info[["UCSC_seqlevel"]][j],
            ":1-",
            chr_info[["UCSC_seqlength"]][j],
            "&sids=", str_c(gtex_metadata_tissue[["rail_id"]], collapse = ",")
          ),
          delim = "\t"
          )
      })
        GTEx_junctions <- bind_rows(GTEx_junctions, junctions_per_chr)
    }

    print(str_c(Sys.time(), " - tidying junction data..."))

    # reformat GTEx junc counts into a count matrix
    GTEx_junctions_tidy <- GTEx_junctions[["samples"]] %>%
      str_split(pattern = ",")

    GTEx_junctions_tidy <- GTEx_junctions_tidy %>%
      lapply(FUN = tidy_raw_count) %>%
      do.call(args = ., bind_rows)

    # format colnames to be R-friendly
    colnames(GTEx_junctions_tidy) <- colnames(GTEx_junctions_tidy) %>% str_c("gtex_", .)

    # replace all NAs with 0s
    GTEx_junctions_tidy[is.na(GTEx_junctions_tidy)] <- 0

    # add back the original junction coordinate details
    GTEx_junctions_tidy <- GTEx_junctions_tidy %>%
      bind_cols(GTEx_junctions %>% dplyr::select(chr = chromosome, start, end, strand), .)

    dest_path <- str_c("data/GTEx_v6_junctions_", tissues_of_interest[["tissue_tidy"]][i], ".rda")

    # upload to dropbox
    save(GTEx_junctions_tidy,
         file = dest_path,
         compress = "gzip"
    )

    drop_upload(dest_path,
                path = "public/dasper/GTEx_v6_junctions"
    )

    file.remove(dest_path)
}

