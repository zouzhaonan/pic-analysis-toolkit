# 役割:
#   biomart 一覧取得、cache 作成、lookup 登録の処理本体をまとめる。
# 入力:
#   source, genome, dataset の CLI 引数。
# 出力:
#   <lib_dir>/register/biomart_lookup.tsv、
#   <lib_dir>/biomart 配下の cache、
#   cache 内 human_ortholog 列。

print_usage <- function() {
  local_help <- file.path(analysis_root, "lib", "help", "help_manage-biomart.txt")
  cat(readChar(local_help, file.info(local_help)$size, useBytes = TRUE))
  cat("\n\n")

  registered <- "  (none)"
  lookup <- tryCatch(read_biomart_lookup(), error = function(err) tibble::tibble())
  if (nrow(lookup) > 0 && "genome" %in% colnames(lookup) && "biomart_dataset" %in% colnames(lookup)) {
    registered_tbl <- lookup |>
      dplyr::filter(!is.na(.data$genome), .data$genome != "", !is.na(.data$biomart_dataset), .data$biomart_dataset != "") |>
      dplyr::distinct(.data$genome, .data$biomart_dataset) |>
      dplyr::arrange(.data$genome)
    if (nrow(registered_tbl) > 0) {
      registered <- paste(sprintf("  - %s (%s)", registered_tbl$genome, registered_tbl$biomart_dataset), collapse = "\n")
    }
  }

  cat("Registered biomart genomes:\n")
  cat(registered)
  cat("\n")
}

parse_args <- function(args) {
  parsed <- list(genome = NULL, dataset_name = NULL, list_datasets = FALSE, register = FALSE, source = "ensembl")
  i <- 1

  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--help")) {
      print_usage()
      quit(save = "no", status = 0)
    }
    if (!startsWith(arg, "--")) {
      stop(sprintf("Unexpected argument: %s", arg), call. = FALSE)
    }

    key <- sub("^--", "", arg)
    if (key %in% c("list-datasets", "register")) {
      parsed[[gsub("-", "_", key)]] <- TRUE
      i <- i + 1
      next
    }

    i <- i + 1
    if (i > length(args)) {
      stop(sprintf("Option requires a value: --%s", key), call. = FALSE)
    }
    value <- args[[i]]

    if (key == "genome") {
      parsed$genome <- value
    } else if (key == "dataset") {
      parsed$dataset_name <- as.character(value)
    } else if (key == "source") {
      parsed$source <- value
    } else {
      stop(sprintf("Unknown option: --%s", key), call. = FALSE)
    }

    i <- i + 1
  }

  if (!parsed$list_datasets && !parsed$register) {
    stop("Specify either --list-datasets or --register.", call. = FALSE)
  }
  if (parsed$list_datasets && (parsed$register || !is.null(parsed$genome) || !is.null(parsed$dataset_name))) {
    stop("--list-datasets cannot be combined with other options.", call. = FALSE)
  }
  if (parsed$register && (is.null(parsed$genome) || is.null(parsed$dataset_name) || parsed$dataset_name == "")) {
    stop("--register requires both --genome and --dataset.", call. = FALSE)
  }
  if (!parsed$register && (!is.null(parsed$genome) || !is.null(parsed$dataset_name))) {
    stop("Use --register together with --genome and --dataset.", call. = FALSE)
  }

  parsed
}

curl_retry_count <- function() {
  4L
}

curl_retry_sleep_sec <- function() {
  3
}

register_retry_count <- function() {
  10L
}

register_retry_sleep_sec <- function() {
  3
}

load_libraries <- function() {
  suppressPackageStartupMessages(library(readr))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tibble))
}

run_curl_capture <- function(args) {
  attempts <- curl_retry_count()
  wait_sec <- curl_retry_sleep_sec()
  last_output <- character()

  for (attempt in seq_len(attempts)) {
    output <- suppressWarnings(system2("curl", args = args, stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    if (is.null(status) || status == 0) {
      return(paste(output, collapse = "\n"))
    }

    last_output <- output
    if (attempt < attempts) {
      message(sprintf("curl capture failed (attempt %d/%d). Retrying in %ss ...", attempt, attempts, wait_sec))
      Sys.sleep(wait_sec)
    }
  }

  stop(paste(last_output, collapse = "\n"), call. = FALSE)
}

run_curl_download <- function(args, output_file) {
  attempts <- curl_retry_count()
  wait_sec <- curl_retry_sleep_sec()
  last_output <- character()

  for (attempt in seq_len(attempts)) {
    output <- suppressWarnings(system2("curl", args = c(args, "-o", output_file), stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    if (is.null(status) || status == 0) {
      return(invisible(NULL))
    }

    last_output <- output
    if (attempt < attempts) {
      message(sprintf("curl download failed (attempt %d/%d). Retrying in %ss ...", attempt, attempts, wait_sec))
      Sys.sleep(wait_sec)
    }
  }

  stop(paste(last_output, collapse = "\n"), call. = FALSE)
}

parse_dataset_table <- function(raw_text, mart_name, host_name) {
  if (!nzchar(trimws(raw_text))) {
    return(tibble::tibble())
  }

  table <- utils::read.delim(text = raw_text, sep = "\t", header = FALSE, quote = "", comment.char = "", fill = TRUE)
  if (nrow(table) == 0) {
    return(tibble::tibble())
  }
  if (ncol(table) == 1 && all(grepl("^Problem retrieving datasets", table[[1]]))) {
    stop(
      sprintf(
        "Failed to retrieve datasets from %s (%s): %s",
        mart_name,
        host_name,
        table[[1]][[1]]
      ),
      call. = FALSE
    )
  }

  while (ncol(table) < 3) {
    table[[ncol(table) + 1]] <- NA_character_
  }

  table <- table[apply(table, 1, function(row) any(nzchar(trimws(as.character(row))))), , drop = FALSE]
  if (nrow(table) == 0) {
    return(tibble::tibble())
  }

  colnames(table)[1:min(9, ncol(table))] <- c(
    "record_type",
    "dataset",
    "description",
    "visible",
    "assembly",
    "date_created",
    "row_limit",
    "default_interface",
    "updated_at"
  )[seq_len(min(9, ncol(table)))]

  tibble::as_tibble(table) |>
    dplyr::filter(.data$record_type == "TableSet") |>
    dplyr::select("dataset", "description", dplyr::everything()) |>
    dplyr::arrange(.data$description, .data$dataset) |>
    dplyr::mutate(biomart = mart_name, host_name = host_name, .before = 1)
}

fetch_available_datasets <- function(source_name) {
  source_info <- resolve_source(source_name)
  mart_name <- source_info$biomart
  host_name <- source_info$host_name

  raw_text <- run_curl_capture(
    c(
      "-L",
      "--silent",
      "--show-error",
      "--get",
      "--data-urlencode",
      "type=datasets",
      "--data-urlencode",
      paste0("mart=", mart_name),
      sprintf("%s/biomart/martservice", host_name)
    )
  )

  parse_dataset_table(raw_text, mart_name, host_name) |>
    dplyr::select("biomart", "host_name", "dataset", "description", dplyr::everything()) |>
    dplyr::mutate(index = dplyr::row_number(), .before = 1)
}

print_available_datasets <- function(source_name) {
  datasets <- fetch_available_datasets(source_name)
  if (nrow(datasets) == 0) {
    stop("No biomart datasets were returned.", call. = FALSE)
  }

  output <- as.data.frame(datasets[, c("index", "dataset", "description")], stringsAsFactors = FALSE)
  utils::write.table(output, file = stdout(), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

normalize_cache_table <- function(cache_file) {
  cache_table <- read_tsv_as_chr(cache_file)
  if (nrow(cache_table) == 0) {
    raw_head <- readLines(cache_file, n = 1, warn = FALSE)
    raw_head <- if (length(raw_head) == 0) "(no content)" else raw_head[[1]]
    stop(
      sprintf("Downloaded biomart cache has no data rows: %s (response head: %s)", cache_file, raw_head),
      call. = FALSE
    )
  }
  if (ncol(cache_table) < 1) {
    stop(sprintf("Downloaded biomart cache does not have required columns: %s", cache_file), call. = FALSE)
  }

  while (ncol(cache_table) < 5) {
    cache_table[[ncol(cache_table) + 1]] <- NA_character_
  }

  cache_table <- cache_table[, 1:5]
  colnames(cache_table) <- c("ens_gene", "ext_gene", "biotype", "chr", "human_ortholog")
  readr::write_tsv(cache_table, cache_file)
}

cache_has_annotation_values <- function(cache_file) {
  cache_table <- tryCatch(read_tsv_as_chr(cache_file), error = function(err) NULL)
  if (is.null(cache_table) || nrow(cache_table) == 0 || ncol(cache_table) < 4) {
    return(FALSE)
  }

  biotype_values <- cache_table[[3]]
  chr_values <- cache_table[[4]]
  has_biotype <- any(!is.na(biotype_values) & biotype_values != "")
  has_chr <- any(!is.na(chr_values) & chr_values != "")
  isTRUE(has_biotype) && isTRUE(has_chr)
}

normalize_symbol_upper <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  x[x %in% c("NA", "N/A")] <- ""
  x[is.na(x)] <- ""
  x
}

attach_human_ortholog_column <- function(cache_file, is_human = FALSE) {
  cache_table <- read_tsv_as_chr(cache_file)
  while (ncol(cache_table) < 5) {
    cache_table[[ncol(cache_table) + 1]] <- ""
  }
  cache_table <- cache_table[, 1:5]
  colnames(cache_table) <- c("ens_gene", "ext_gene", "biotype", "chr", "human_ortholog")
  cache_table <- cache_table |>
    dplyr::mutate(
      ens_gene = as.character(.data$ens_gene),
      ext_gene = as.character(.data$ext_gene),
      biotype = as.character(.data$biotype),
      chr = as.character(.data$chr),
      human_ortholog = as.character(.data$human_ortholog)
    )

  if (isTRUE(is_human)) {
    cache_table <- cache_table |>
      dplyr::mutate(human_ortholog = normalize_symbol_upper(.data$ext_gene))
    readr::write_tsv(cache_table, cache_file)
    return(invisible(NULL))
  }

  cache_table <- cache_table |>
    dplyr::mutate(human_ortholog = normalize_symbol_upper(.data$human_ortholog)) |>
    dplyr::select("ens_gene", "ext_gene", "biotype", "chr", "human_ortholog")

  readr::write_tsv(cache_table, cache_file)
}

lookup_dataset_by_name <- function(dataset_name, source_name) {
  datasets <- fetch_available_datasets(source_name)
  selected <- datasets |>
    dplyr::filter(.data$dataset == dataset_name)
  if (nrow(selected) == 0) {
    stop(sprintf("No biomart dataset was found for dataset '%s'.", dataset_name), call. = FALSE)
  }

  selected
}

download_cache_with_host <- function(host_name, dataset_name, cache_file, require_human_ortholog = FALSE) {
  dir.create(biomart_cache_dir, recursive = TRUE, showWarnings = FALSE)
  temp_cache_file <- tempfile(pattern = paste0(dataset_name, "_"), fileext = ".tsv", tmpdir = biomart_cache_dir)
  on.exit(if (file.exists(temp_cache_file)) unlink(temp_cache_file), add = TRUE)
  attribute_sets <- list(
    c("ensembl_gene_id", "external_gene_name", "gene_biotype", "chromosome_name"),
    c("ensembl_gene_id", "external_gene_name", "gene_type", "chromosome_name"),
    c("ensembl_gene_id", "external_gene_name", "gene_biotype", "chromosome_scaffold_name"),
    c("ensembl_gene_id", "external_gene_name", "gene_type", "chromosome_scaffold_name")
  )
  last_error <- NULL
  succeeded <- FALSE
  for (attrs in attribute_sets) {
    query_xml <- paste0(
      '<?xml version="1.0" encoding="UTF-8"?>',
      '<!DOCTYPE Query>',
      '<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6">',
      sprintf('<Dataset name="%s" interface="default">', dataset_name),
      paste(sprintf('<Attribute name="%s"/>', attrs), collapse = ""),
      '</Dataset>',
      '</Query>'
    )

    query_file <- tempfile(pattern = "biomart_query_", fileext = ".xml")
    writeLines(query_xml, query_file, useBytes = TRUE)
    on.exit(if (file.exists(query_file)) unlink(query_file), add = TRUE)

    attempt_ok <- tryCatch(
      {
        run_curl_download(
          c(
            "-L",
            "--silent",
            "--show-error",
            "--get",
            "--data-urlencode",
            paste0("query@", query_file),
            sprintf("%s/biomart/martservice", host_name)
          ),
          temp_cache_file
        )

        if (!file.exists(temp_cache_file) || file.info(temp_cache_file)$size == 0) {
          stop(sprintf("Downloaded biomart cache is empty: %s", temp_cache_file), call. = FALSE)
        }

        normalize_cache_table(temp_cache_file)
        if (!cache_has_annotation_values(temp_cache_file)) {
          stop(
            sprintf(
              "Downloaded biomart cache has no usable biotype/chr values for attrs: %s",
              paste(attrs, collapse = ",")
            ),
            call. = FALSE
          )
        }
        TRUE
      },
      error = function(err) {
        last_error <<- conditionMessage(err)
        FALSE
      }
    )

    if (isTRUE(attempt_ok)) {
      succeeded <- TRUE
      break
    }
  }

  if (!isTRUE(succeeded)) {
    stop(last_error, call. = FALSE)
  }

  if (!file.rename(temp_cache_file, cache_file)) {
    stop(sprintf("Failed to move biomart cache into place: %s", cache_file), call. = FALSE)
  }

  message(sprintf("Saved biomart cache: %s", cache_file))
}

fill_human_ortholog_with_host <- function(host_name, dataset_name, cache_file, rounds = 8L, sleep_sec = 5) {
  if (!file.exists(cache_file)) return(invisible(NULL))

  cache_tbl <- read_tsv_as_chr(cache_file)
  while (ncol(cache_tbl) < 5) cache_tbl[[ncol(cache_tbl) + 1]] <- ""
  cache_tbl <- cache_tbl[, 1:5]
  colnames(cache_tbl) <- c("ens_gene", "ext_gene", "biotype", "chr", "human_ortholog")

  for (round_idx in seq_len(as.integer(rounds))) {
    query_xml <- paste0(
      '<?xml version="1.0" encoding="UTF-8"?>',
      '<!DOCTYPE Query>',
      '<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">',
      sprintf('<Dataset name="%s" interface="default">', dataset_name),
      '<Attribute name="ensembl_gene_id"/>',
      '<Attribute name="hsapiens_homolog_associated_gene_name"/>',
      '</Dataset>',
      '</Query>'
    )
    query_file <- tempfile(pattern = "biomart_hs_ortholog_query_", fileext = ".xml")
    writeLines(query_xml, query_file, useBytes = TRUE)
    on.exit(if (file.exists(query_file)) unlink(query_file), add = TRUE)

    tmp_file <- tempfile(pattern = "biomart_hs_ortholog_", fileext = ".tsv")
    on.exit(if (file.exists(tmp_file)) unlink(tmp_file), add = TRUE)

    ok <- tryCatch(
      {
        run_curl_download(
          c(
            "-L",
            "--silent",
            "--show-error",
            "--get",
            "--data-urlencode",
            paste0("query@", query_file),
            sprintf("%s/biomart/martservice", host_name)
          ),
          tmp_file
        )
        TRUE
      },
      error = function(err) FALSE
    )

    if (isTRUE(ok) && file.exists(tmp_file) && file.info(tmp_file)$size > 0) {
      orth_tbl <- read_tsv_as_chr(tmp_file)
      if (ncol(orth_tbl) >= 2) {
        orth_tbl <- orth_tbl[, 1:2]
        colnames(orth_tbl) <- c("ens_gene", "human_ortholog")
        orth_tbl <- orth_tbl |>
          dplyr::transmute(
            ens_gene = as.character(.data$ens_gene),
            human_ortholog = normalize_symbol_upper(.data$human_ortholog)
          ) |>
          dplyr::filter(!is.na(.data$ens_gene), .data$ens_gene != "", !is.na(.data$human_ortholog), .data$human_ortholog != "") |>
          dplyr::distinct(.data$ens_gene, .data$human_ortholog)

        if (nrow(orth_tbl) > 0) {
          cache_tbl <- cache_tbl |>
            dplyr::mutate(
              ens_gene = as.character(.data$ens_gene),
              human_ortholog = normalize_symbol_upper(.data$human_ortholog)
            ) |>
            dplyr::left_join(orth_tbl, by = "ens_gene") |>
            dplyr::mutate(
              human_ortholog = dplyr::if_else(
                !is.na(.data$human_ortholog.y) & .data$human_ortholog.y != "",
                .data$human_ortholog.y,
                .data$human_ortholog.x
              )
            ) |>
            dplyr::select("ens_gene", "ext_gene", "biotype", "chr", "human_ortholog")
          readr::write_tsv(cache_tbl, cache_file)
          mapped_n <- cache_tbl |>
            dplyr::summarise(n = sum(!is.na(.data$human_ortholog) & .data$human_ortholog != "" & .data$human_ortholog != "NA")) |>
            dplyr::pull(.data$n)
          if (isTRUE(mapped_n > 0)) {
            message(sprintf("human_ortholog mapped rows: %d", as.integer(mapped_n)))
            return(invisible(NULL))
          }
        }
      }
    }

    if (round_idx < rounds) {
      message(sprintf("human_ortholog not available yet (round %d/%d). Retrying in %ss ...", round_idx, rounds, sleep_sec))
      Sys.sleep(sleep_sec)
    }
  }

  stop(sprintf("Failed to fetch human_ortholog values after %d rounds: %s", rounds, cache_file), call. = FALSE)
}

register_lookup_row <- function(genome, biomart_name, dataset_name, description, host_name) {
  current_lookup <- read_biomart_lookup()

  new_row <- tibble::tibble(
    genome = genome,
    biomart = biomart_name,
    biomart_dataset = dataset_name,
    description = description,
    host_name = host_name
  )

  replaced <- any(current_lookup$genome == genome)
  if (replaced) {
    current_lookup <- current_lookup |>
      dplyr::filter(.data$genome != genome)
  }

  updated_lookup <- dplyr::bind_rows(current_lookup, new_row)
  readr::write_tsv(updated_lookup, biomart_lookup_file)
  if (replaced) {
    message(sprintf("Updated genome lookup: %s -> %s", genome, dataset_name))
  } else {
    message(sprintf("Registered genome lookup: %s -> %s", genome, dataset_name))
  }
}

is_biomart_registration_complete <- function(genome, dataset_name, cache_file) {
  cache_ok <- file.exists(cache_file) && file.info(cache_file)$size > 0
  if (!isTRUE(cache_ok)) {
    return(FALSE)
  }

  lookup_ok <- FALSE
  lookup <- tryCatch(read_biomart_lookup(), error = function(err) tibble::tibble())
  if (nrow(lookup) > 0 && all(c("genome", "biomart_dataset") %in% colnames(lookup))) {
    lookup_ok <- any(lookup$genome == genome & lookup$biomart_dataset == dataset_name)
  }
  isTRUE(lookup_ok)
}

register_biomart_with_retry <- function(genome, dataset_name, source_name) {
  selected_dataset <- lookup_dataset_by_name(dataset_name, source_name)
  dataset_name <- selected_dataset$dataset[[1]]
  description <- selected_dataset$description[[1]]
  biomart_name <- selected_dataset$biomart[[1]]
  host_name <- selected_dataset$host_name[[1]]
  cache_file <- file.path(biomart_cache_dir, paste0(dataset_name, ".tsv"))

  attempts <- register_retry_count()
  wait_sec <- register_retry_sleep_sec()
  last_error <- NULL

  for (attempt in seq_len(attempts)) {
    tryCatch(
      {
        if (file.exists(cache_file)) {
          message(sprintf("Overwriting existing biomart cache: %s", cache_file))
        }

        download_cache_with_host(host_name, dataset_name, cache_file, require_human_ortholog = FALSE)
        if (dataset_name != "hsapiens_gene_ensembl") {
          fill_human_ortholog_with_host(host_name, dataset_name, cache_file)
          attach_human_ortholog_column(cache_file, is_human = FALSE)
        } else {
          attach_human_ortholog_column(cache_file, is_human = TRUE)
        }

        register_lookup_row(genome, biomart_name, dataset_name, description, host_name)
      },
      error = function(err) {
        last_error <<- conditionMessage(err)
      }
    )

    if (is_biomart_registration_complete(genome, dataset_name, cache_file)) {
      message(sprintf("biomart registration completed: %s -> %s", genome, dataset_name))
      return(invisible(NULL))
    }

    if (attempt < attempts) {
      message(sprintf("biomart registration failed (attempt %d/%d). Retrying in %ss ...", attempt, attempts, wait_sec))
      Sys.sleep(wait_sec)
    }
  }

  if (is.null(last_error) || !nzchar(last_error)) {
    last_error <- sprintf("cache/lookup verification failed for genome=%s dataset=%s", genome, dataset_name)
  }
  stop(
    sprintf("Failed to register biomart after %d attempts: %s", attempts, last_error),
    call. = FALSE
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  load_libraries()

  if (args$list_datasets) {
    print_available_datasets(args$source)
    return(invisible(NULL))
  }

  if (args$register) {
    register_biomart_with_retry(
      genome = args$genome,
      dataset_name = args$dataset_name,
      source_name = args$source
    )
    return(invisible(NULL))
  }
}
