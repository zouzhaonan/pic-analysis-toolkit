# 役割:
#   DESeq2 stats から enrichment を実行する補助関数を提供する。
#   enrichment 資材は local gson を利用し、非対応種は human gson へフォールバックする。

normalize_method_names <- function(methods) {
  x <- toupper(trimws(as.character(methods)))
  x <- dplyr::recode(x, MGI_MPO = "MPO", .default = x)
  x
}

supported_enrichment_methods <- function() {
  c("GO_BP", "GO_CC", "GO_MF", "HPO", "HDO", "MPO", "KEGG", "REACTOME", "WIKIPATHWAYS")
}

validate_ora_methods <- function(methods) {
  methods <- normalize_method_names(methods)
  unsupported <- setdiff(methods, supported_enrichment_methods())
  if (length(unsupported) > 0) {
    stop(
      sprintf(
        "Unsupported ORA method(s): %s. Supported methods: %s",
        paste(unsupported, collapse = ", "),
        paste(supported_enrichment_methods(), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  methods
}

validate_gsea_methods <- function(methods) {
  methods <- normalize_method_names(methods)
  unsupported <- setdiff(methods, supported_enrichment_methods())
  if (length(unsupported) > 0) {
    stop(
      sprintf(
        "Unsupported GSEA method(s): %s. Supported methods: %s",
        paste(unsupported, collapse = ", "),
        paste(supported_enrichment_methods(), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  methods
}

has_text_value <- function(x) {
  !is.null(x) && length(x) > 0 && !is.na(x[[1]]) && nzchar(x[[1]])
}

parse_bool_flag <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x %in% c("true", "t", "1", "yes", "y")
}

normalize_human_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  x[x %in% c("", "NA", "N/A")] <- NA_character_
  x
}

normalize_source_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A")] <- NA_character_
  x
}

safe_as_tibble <- function(result_obj) {
  if (is.null(result_obj)) {
    return(tibble::tibble())
  }
  tibble::as_tibble(as.data.frame(result_obj))
}

read_genome_enrichment_map <- function(genome_map_file) {
  if (!has_text_value(genome_map_file) || !file.exists(genome_map_file)) {
    return(tibble::tibble())
  }
  suppressMessages(
    readr::read_tsv(
      genome_map_file,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE
    )
  )
}

resolve_gson_row <- function(genome, method, genome_map_file) {
  method <- toupper(as.character(method))
  map_tbl <- read_genome_enrichment_map(genome_map_file)
  if (nrow(map_tbl) == 0) return(NULL)

  # Backward compatibility: old long format.
  if (all(c("genome", "resource", "species_key") %in% colnames(map_tbl))) {
    own <- map_tbl |>
      dplyr::filter(.data$genome == .env$genome, toupper(.data$resource) == .env$method) |>
      dplyr::slice_head(n = 1)
    if (nrow(own) > 0) return(own)
    human <- map_tbl |>
      dplyr::filter(.data$genome == "hg38", toupper(.data$resource) == .env$method) |>
      dplyr::slice_head(n = 1)
    if (nrow(human) > 0) return(human)
    return(NULL)
  }

  # New wide format: genome x category boolean.
  if (!all(c("genome", "species_key") %in% colnames(map_tbl))) return(NULL)
  if (!method %in% colnames(map_tbl)) return(NULL)

  own <- map_tbl |>
    dplyr::filter(.data$genome == .env$genome) |>
    dplyr::slice_head(n = 1)

  if (nrow(own) > 0 && has_text_value(as.character(own$species_key[[1]])) && parse_bool_flag(own[[method]][[1]])) {
    return(tibble::tibble(genome = as.character(own$genome[[1]]), species_key = as.character(own$species_key[[1]])))
  }

  # Fallback to human gson when species-specific resource is unavailable.
  human <- map_tbl |>
    dplyr::filter(.data$genome == "hg38") |>
    dplyr::slice_head(n = 1)
  if (nrow(human) > 0 && has_text_value(as.character(human$species_key[[1]])) && parse_bool_flag(human[[method]][[1]])) {
    return(tibble::tibble(genome = as.character(human$genome[[1]]), species_key = as.character(human$species_key[[1]])))
  }

  NULL
}

is_method_available_for_genome <- function(genome, method, genome_map_file) {
  method <- toupper(as.character(method))
  map_tbl <- read_genome_enrichment_map(genome_map_file)
  if (nrow(map_tbl) == 0) return(FALSE)

  # Backward compatibility: old long format.
  if (all(c("genome", "resource", "species_key") %in% colnames(map_tbl))) {
    own <- map_tbl |>
      dplyr::filter(.data$genome == .env$genome, toupper(.data$resource) == .env$method) |>
      dplyr::slice_head(n = 1)
    return(nrow(own) > 0)
  }

  if (!all(c("genome", "species_key") %in% colnames(map_tbl))) return(FALSE)
  if (!method %in% colnames(map_tbl)) return(FALSE)

  own <- map_tbl |>
    dplyr::filter(.data$genome == .env$genome) |>
    dplyr::slice_head(n = 1)
  if (nrow(own) == 0) return(FALSE)
  isTRUE(parse_bool_flag(own[[method]][[1]]))
}

resolve_enrichment_gene_column <- function(genome, method, genome_map_file) {
  if (isTRUE(is_method_available_for_genome(genome, method, genome_map_file))) {
    return("source_gene_symbol")
  }
  "human_gene_symbol"
}

read_local_gson <- function(genome, method, genome_map_file, enrichment_lib_dir) {
  row <- resolve_gson_row(genome = genome, method = method, genome_map_file = genome_map_file)
  if (is.null(row) || nrow(row) == 0) return(NULL)

  species_key <- as.character(row$species_key[[1]])
  if (!has_text_value(species_key)) return(NULL)
  path <- file.path(enrichment_lib_dir, sprintf("%s__%s.gson", toupper(as.character(method)), species_key))
  if (!has_text_value(path) || !file.exists(path)) return(NULL)
  readRDS(path)
}

run_local_ora <- function(genes_symbol, universe_symbol, gson_obj) {
  if (is.null(gson_obj) || length(genes_symbol) == 0) return(NULL)
  quiet_cluster_call <- function(f) {
    out <- NULL
    suppressWarnings(
      suppressMessages(
        utils::capture.output(
          utils::capture.output(
            {
              out <- f()
            },
            type = "message"
          ),
          type = "output"
        )
      )
    )
    out
  }
  tryCatch(
    quiet_cluster_call(function() {
      clusterProfiler::enricher(
        gene = genes_symbol,
        universe = universe_symbol,
        gson = gson_obj,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )
    }),
    error = function(err) NULL
  )
}

run_local_gsea <- function(gene_list_symbol, gson_obj) {
  if (is.null(gson_obj) || length(gene_list_symbol) == 0) return(NULL)
  quiet_cluster_call <- function(f) {
    out <- NULL
    suppressWarnings(
      suppressMessages(
        utils::capture.output(
          utils::capture.output(
            {
              out <- f()
            },
            type = "message"
          ),
          type = "output"
        )
      )
    )
    out
  }
  tryCatch(
    quiet_cluster_call(function() {
      clusterProfiler::GSEA(
        geneList = gene_list_symbol,
        gson = gson_obj,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 1,
        verbose = FALSE
      )
    }),
    error = function(err) NULL
  )
}

read_tsv_as_chr_local <- function(path) {
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

resolve_biomart_lookup_row <- function(genome, biomart_lookup_file) {
  lookup <- read_tsv_as_chr_local(biomart_lookup_file) |>
    dplyr::filter(.data$genome == .env$genome)
  if (nrow(lookup) == 0) {
    stop(
      sprintf(
        "No biomart lookup row for genome '%s'. Run: pic manage-biomart --register --genome %s --dataset <dataset_name>",
        genome, genome
      ),
      call. = FALSE
    )
  }
  lookup |> dplyr::slice_head(n = 1)
}

load_source_to_human_symbol_map <- function(genome, biomart_lookup_file, biomart_cache_dir, biomart_ortholog_dir = NULL) {
  row <- resolve_biomart_lookup_row(genome, biomart_lookup_file)
  dataset <- as.character(row$biomart_dataset[[1]])
  cache_file <- file.path(biomart_cache_dir, paste0(dataset, ".tsv"))
  if (!file.exists(cache_file)) {
    stop(sprintf("biomart cache was not found: %s", cache_file), call. = FALSE)
  }

  cache_df <- read_tsv_as_chr_local(cache_file)
  if ("human_ortholog" %in% colnames(cache_df)) {
    from_cache <- cache_df |>
      dplyr::transmute(
        source_gene_id = as.character(.data$ens_gene),
        human_gene_symbol = normalize_human_symbol(.data$human_ortholog)
      ) |>
      dplyr::filter(!is.na(.data$source_gene_id), .data$source_gene_id != "", !is.na(.data$human_gene_symbol), .data$human_gene_symbol != "") |>
      dplyr::distinct(.data$source_gene_id, .data$human_gene_symbol)
    if (nrow(from_cache) > 0) {
      return(from_cache)
    }
  }

  if (!is.na(dataset) && dataset == "hsapiens_gene_ensembl") {
    map_df <- cache_df
    if (!all(c("ens_gene", "ext_gene") %in% colnames(map_df))) {
      colnames(map_df)[1:min(2, ncol(map_df))] <- c("ens_gene", "ext_gene")[seq_len(min(2, ncol(map_df)))]
    }
    return(
      map_df |>
        dplyr::transmute(
          source_gene_id = as.character(.data$ens_gene),
          human_gene_symbol = normalize_human_symbol(.data$ext_gene)
        ) |>
        dplyr::filter(!is.na(.data$source_gene_id), .data$source_gene_id != "", !is.na(.data$human_gene_symbol), .data$human_gene_symbol != "") |>
        dplyr::distinct(.data$source_gene_id, .data$human_gene_symbol)
    )
  }

  cache_df |>
    dplyr::transmute(
      source_gene_id = as.character(.data$ens_gene),
      human_gene_symbol = normalize_human_symbol(.data$ext_gene)
    ) |>
    dplyr::filter(!is.na(.data$source_gene_id), .data$source_gene_id != "", !is.na(.data$human_gene_symbol), .data$human_gene_symbol != "") |>
    dplyr::distinct(.data$source_gene_id, .data$human_gene_symbol)
}

map_stats_to_human_symbols <- function(stats, genome, biomart_lookup_file, biomart_cache_dir, biomart_ortholog_dir) {
  if (!"ens_gene" %in% colnames(stats)) {
    stop("stats table must include 'ens_gene' column.", call. = FALSE)
  }

  id_map <- load_source_to_human_symbol_map(genome, biomart_lookup_file, biomart_cache_dir, biomart_ortholog_dir) |>
    dplyr::filter(!is.na(.data$source_gene_id), .data$source_gene_id != "", !is.na(.data$human_gene_symbol), .data$human_gene_symbol != "") |>
    dplyr::arrange(.data$source_gene_id, .data$human_gene_symbol) |>
    dplyr::distinct(.data$source_gene_id, .data$human_gene_symbol)

  mapped <- stats |>
    dplyr::left_join(
      id_map,
      by = c("ens_gene" = "source_gene_id"),
      relationship = "many-to-many"
    )

  mapped <- mapped |>
    dplyr::mutate(
      ens_gene = as.character(.data$ens_gene),
      source_gene_symbol = if ("ext_gene" %in% colnames(mapped)) normalize_source_symbol(.data$ext_gene) else NA_character_,
      human_gene_symbol = normalize_human_symbol(.data$human_gene_symbol)
    )

  if ("ext_gene" %in% colnames(mapped)) {
    mapped <- mapped |>
      dplyr::mutate(
        human_gene_symbol = dplyr::if_else(
          !is.na(.data$human_gene_symbol) & .data$human_gene_symbol != "",
          .data$human_gene_symbol,
          normalize_human_symbol(.data$ext_gene)
        )
      )
  }

  mapped |>
    dplyr::mutate(
      human_gene_symbol = dplyr::if_else(
        !is.na(.data$human_gene_symbol) & .data$human_gene_symbol != "",
        .data$human_gene_symbol,
        normalize_human_symbol(.data$ens_gene)
      ),
      source_gene_symbol = dplyr::if_else(
        !is.na(.data$source_gene_symbol) & .data$source_gene_symbol != "",
        .data$source_gene_symbol,
        normalize_source_symbol(.data$ens_gene)
      )
    )
}

run_single_ora_enrichment <- function(method, genes, universe, genome, genome_map_file, enrichment_lib_dir = NULL) {
  gson_obj <- read_local_gson(genome = genome, method = method, genome_map_file = genome_map_file, enrichment_lib_dir = enrichment_lib_dir)
  res <- run_local_ora(genes, universe, gson_obj)
  safe_as_tibble(res)
}

run_single_gsea_enrichment <- function(method, gene_list, genome, genome_map_file, enrichment_lib_dir = NULL) {
  gson_obj <- read_local_gson(genome = genome, method = method, genome_map_file = genome_map_file, enrichment_lib_dir = enrichment_lib_dir)
  res <- run_local_gsea(gene_list, gson_obj)
  safe_as_tibble(res)
}
