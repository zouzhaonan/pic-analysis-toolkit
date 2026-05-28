# 役割:
#   DESeq2 側で biomart lookup と cache を読む処理をまとめる。
# 入力:
#   genome 名、biomart_lookup_file。
# 出力:
#   ens_gene, ext_gene, biotype, chr, human_ortholog の tibble。

source(file.path(script_dir, "biomart_shared.R"))

load_e2g_from_cache <- function(genome_name, biomart_lookup_file) {
  lookup <- read_tsv_as_chr(biomart_lookup_file) |>
    dplyr::filter(.data$genome == genome_name)

  if (nrow(lookup) == 0 || is.na(lookup$biomart_dataset[[1]]) || lookup$biomart_dataset[[1]] == "") {
    stop(
      sprintf(
        paste(
          "No biomart setting was found for genome '%s'.",
          "Run `pic manage-biomart --list-datasets`, choose a dataset,",
          "then run `pic manage-biomart --register --genome %s --dataset <dataset_name>`.",
          "After that, run `pic deseq2` again.",
          sep = "\n"
        ),
        genome_name,
        genome_name
      ),
      call. = FALSE
    )
  }

  cache_file <- file.path(
    biomart_cache_dir,
    paste0(as.character(lookup$biomart_dataset[[1]]), ".tsv")
  )
  if (!file.exists(cache_file)) {
    stop(
      sprintf(
        paste(
          "biomart cache was not found: %s",
          "Run `pic manage-biomart --list-datasets`, choose a dataset,",
          "then run `pic manage-biomart --register --genome %s --dataset <dataset_name>`.",
          "After that, run `pic deseq2` again.",
          sep = "\n"
        ),
        cache_file,
        genome_name
      ),
      call. = FALSE
    )
  }

  e2g <- read_tsv_as_chr(cache_file)
  while (ncol(e2g) < 5) {
    e2g[[ncol(e2g) + 1]] <- NA_character_
  }
  e2g <- e2g[, 1:5]
  colnames(e2g) <- c("ens_gene", "ext_gene", "biotype", "chr", "human_ortholog")
  tibble::as_tibble(e2g)
}
