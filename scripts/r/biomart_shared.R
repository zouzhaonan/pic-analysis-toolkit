# 役割:
#   biomart lookup と cache の共通設定と補助関数をまとめる。
# 入力:
#   呼び出し元で設定した biomart_lookup_file。
# 出力:
#   lookup 読み込み、source 解決、cache 参照に使う関数。

biomart_cache_dir <- pic_runtime_paths(analysis_root)$biomart_cache_dir
biomart_ortholog_dir <- pic_runtime_paths(analysis_root)$biomart_ortholog_dir

biomart_sources <- list(
  ensembl = list(biomart = "ENSEMBL_MART_ENSEMBL", host_name = "https://www.ensembl.org"),
  plants = list(biomart = "plants_mart", host_name = "https://plants.ensembl.org"),
  fungi = list(biomart = "fungi_mart", host_name = "https://fungi.ensembl.org"),
  protists = list(biomart = "protists_mart", host_name = "https://protists.ensembl.org"),
  metazoa = list(biomart = "metazoa_mart", host_name = "https://metazoa.ensembl.org")
)

read_tsv_as_chr <- function(path) {
  tbl <- suppressMessages(
    readr::read_tsv(
      path,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE
    )
  )
  for (nm in colnames(tbl)) {
    tbl[[nm]] <- gsub("\r", "", as.character(tbl[[nm]]), fixed = TRUE)
  }
  tbl
}

read_biomart_lookup <- function() {
  if (!file.exists(biomart_lookup_file)) {
    dir.create(dirname(biomart_lookup_file), recursive = TRUE, showWarnings = FALSE)
    empty_lookup <- tibble::tibble(
      genome = character(),
      biomart = character(),
      biomart_dataset = character(),
      description = character(),
      host_name = character()
    )
    readr::write_tsv(empty_lookup, biomart_lookup_file)
  }

  tbl <- read_tsv_as_chr(biomart_lookup_file)
  required_cols <- c("genome", "biomart", "biomart_dataset", "description", "host_name")
  missing_cols <- setdiff(required_cols, colnames(tbl))
  if (length(missing_cols) > 0) {
    for (col_name in missing_cols) {
      tbl[[col_name]] <- ""
    }
    tbl <- tbl |>
      dplyr::select(dplyr::all_of(required_cols), dplyr::everything())
    readr::write_tsv(tbl, biomart_lookup_file)
  }

  tbl
}

resolve_source <- function(source_name) {
  source_key <- tolower(source_name)
  source_info <- biomart_sources[[source_key]]
  if (is.null(source_info)) {
    stop(
      sprintf(
        "Unknown source: %s\nAvailable sources: %s",
        source_name,
        paste(names(biomart_sources), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  source_info
}
