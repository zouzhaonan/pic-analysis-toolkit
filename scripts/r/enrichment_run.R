run_jobs_parallel <- function(jobs, parallel_threads, fn) {
  if (length(jobs) == 0) {
    return(invisible(NULL))
  }

  threads <- suppressWarnings(as.integer(parallel_threads))
  if (is.na(threads) || threads < 1) {
    threads <- 1L
  }

  worker <- function(i) {
    tryCatch(
      {
        fn(jobs[[i]])
        list(ok = TRUE, idx = i, error = NULL)
      },
      error = function(err) list(ok = FALSE, idx = i, error = conditionMessage(err))
    )
  }

  results <- if (.Platform$OS.type == "unix" && threads > 1) {
    parallel::mclapply(
      seq_along(jobs),
      worker,
      mc.cores = min(threads, length(jobs))
    )
  } else {
    lapply(seq_along(jobs), worker)
  }

  failed <- Filter(function(x) !isTRUE(x$ok), results)
  if (length(failed) > 0) {
    details <- vapply(
      failed,
      function(x) {
        job <- jobs[[x$idx]]
        label <- if (!is.null(job$contrast_name) && length(job$contrast_name) > 0) {
          as.character(job$contrast_name)
        } else if (!is.null(job$cluster_id) && length(job$cluster_id) > 0) {
          as.character(job$cluster_id)
        } else {
          sprintf("job_%d", as.integer(x$idx))
        }
        sprintf(
          "method=%s contrast=%s error=%s",
          as.character(job$method),
          label,
          as.character(x$error)
        )
      },
      character(1)
    )
    stop(
      sprintf(
        "Enrichment parallel jobs failed (%d/%d). %s",
        length(failed),
        length(jobs),
        paste(head(details, 5), collapse = " | ")
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

strip_kegg_species_suffix <- function(description_vec) {
  x <- as.character(description_vec)
  x <- trimws(x)
  # e.g. "Insulin secretion - Mus musculus (house mouse)" -> "Insulin secretion"
  sub("\\s+-\\s+[^-]+\\([^)]*\\)\\s*$", "", x, perl = TRUE)
}

method_output_meta <- function(method, mode = c("ora", "gsea")) {
  mode <- match.arg(mode)

  if (mode == "ora") {
    return(switch(
      method,
      GO_BP = list(prefix = "GO_BP"),
      GO_CC = list(prefix = "GO_CC"),
      GO_MF = list(prefix = "GO_MF"),
      list(prefix = method)
    ))
  }

  switch(
    method,
    GO_BP = list(prefix = "GSEA_GO_BP"),
    GO_CC = list(prefix = "GSEA_GO_CC"),
    GO_MF = list(prefix = "GSEA_GO_MF"),
    list(prefix = paste0("GSEA_", method))
  )
}

# enrich 出力は CSV のみ (静的画像は廃止)。csv/ ラッパを外し
# out_dir/ORA/<method>/, out_dir/GSEA/<method>/ に直接置く。
ora_csv_dir <- function(out_dir, method) {
  file.path(out_dir, "ORA", method)
}

gsea_csv_dir <- function(out_dir, method) {
  file.path(out_dir, "GSEA", method)
}


run_ora_enrichment <- function(stats, genome, out_dir, project_name,
                               methods = c("GO_BP"), parallel_threads = 1,
                               enrichment_lib_dir = NULL,
                               genome_map_file = NULL,
                               deg_cluster_file = NULL) {
  methods <- validate_ora_methods(methods)
  if (length(methods) == 0) {
    return(invisible(NULL))
  }
  if (is.null(deg_cluster_file) || !nzchar(as.character(deg_cluster_file))) {
    stop("ORA requires --deg-clusters.", call. = FALSE)
  }

  purrr::walk(methods, function(method) {
    method_dir <- ora_csv_dir(out_dir, method)
    unlink(method_dir, recursive = TRUE)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
  })

  build_method_clusters <- function(method_gene_col, method_precomputed_clusters = NULL) {
    method_precomputed_clusters
  }

  first_method_gene_col <- resolve_enrichment_gene_column(
    genome = genome,
    method = methods[[1]],
    genome_map_file = genome_map_file
  )
  first_precomputed_clusters <- load_precomputed_deg_clusters(
    cluster_file = deg_cluster_file,
    stats = stats,
    gene_col = first_method_gene_col
  )
  first_clusters <- build_method_clusters(
    method_gene_col = first_method_gene_col,
    method_precomputed_clusters = first_precomputed_clusters
  )
  message(sprintf("[INFO] ORA clustering by degPatterns: %d clusters", nrow(first_clusters)))

  purrr::walk(methods, function(method) {
    method_gene_col <- resolve_enrichment_gene_column(genome = genome, method = method, genome_map_file = genome_map_file)
    message(sprintf("[INFO] ORA %s uses gene column: %s", method, method_gene_col))

    method_stats <- stats |>
      dplyr::mutate(enrich_gene_symbol = as.character(.data[[method_gene_col]])) |>
      dplyr::filter(!is.na(.data$enrich_gene_symbol), .data$enrich_gene_symbol != "")

    method_precomputed_clusters <- load_precomputed_deg_clusters(
      cluster_file = deg_cluster_file,
      stats = stats,
      gene_col = method_gene_col
    )

    clusters <- build_method_clusters(method_gene_col, method_precomputed_clusters = method_precomputed_clusters)
    if (nrow(clusters) == 0) {
      message(sprintf("[INFO] ORA %s skipped: no clusters after mapping", method))
      return(invisible(NULL))
    }

    jobs <- list()
    for (idx in seq_len(nrow(clusters))) {
      jobs[[length(jobs) + 1]] <- list(
        method = method,
        cluster_id = as.character(clusters$cluster_id[[idx]]),
        cluster_genes = as.character(clusters$genes[[idx]])
      )
    }

    run_jobs_parallel(jobs, parallel_threads, function(job) {
      cluster_id <- as.character(job$cluster_id)
      cluster_genes <- unique(as.character(job$cluster_genes))
      display_label <- cluster_id
      file_label <- format_contrast_file_label(cluster_id)
      meta <- method_output_meta(method, mode = "ora")
      method_dir <- ora_csv_dir(out_dir, method)

      if (length(cluster_genes) == 0) {
        empty_cols <- c(
          "method", "cluster", "ID", "Description", "GeneRatio", "BgRatio",
          "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count"
        )
        empty_tbl <- as.data.frame(setNames(rep(list(character()), length(empty_cols)), empty_cols), stringsAsFactors = FALSE)
        readr::write_csv(
          tibble::as_tibble(empty_tbl),
          file.path(method_dir, sprintf("%s_%s_%s.csv", meta$prefix, file_label, project_name))
        )
        return(invisible(NULL))
      }

      universe_genes <- method_stats |>
        dplyr::pull(.data$enrich_gene_symbol) |>
        unique()

      cluster_res <- run_single_ora_enrichment(
        method = method,
        genes = cluster_genes,
        universe = universe_genes,
        genome = genome,
        genome_map_file = genome_map_file,
        enrichment_lib_dir = enrichment_lib_dir
      ) |>
        dplyr::mutate(direction = cluster_id)

      required_ora_cols <- c(
        "ID", "Description", "GeneRatio", "BgRatio",
        "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "direction"
      )
      numeric_ora_cols <- c("RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "Count")
      normalize_ora_result <- function(tbl) {
        missing <- setdiff(required_ora_cols, colnames(tbl))
        if (length(missing) > 0) {
          for (col_name in missing) {
            if (col_name %in% numeric_ora_cols) {
              tbl[[col_name]] <- NA_real_
            } else {
              tbl[[col_name]] <- NA_character_
            }
          }
        }
        for (col_name in required_ora_cols) {
          if (col_name %in% numeric_ora_cols) {
            tbl[[col_name]] <- suppressWarnings(as.numeric(tbl[[col_name]]))
          } else {
            tbl[[col_name]] <- as.character(tbl[[col_name]])
          }
        }
        tbl
      }
      cluster_res <- normalize_ora_result(cluster_res)

      result_cols <- c(
        "cluster", "ID", "Description", "GeneRatio", "BgRatio",
        "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count"
      )

      is_kegg_method <- identical(as.character(method), "KEGG")
      result_table <- cluster_res |>
        dplyr::filter(!is.na(.data$pvalue)) |>
        dplyr::arrange(.data$pvalue) |>
        dplyr::mutate(
          cluster = display_label,
          method = method,
          Description = if (is_kegg_method) strip_kegg_species_suffix(.data$Description) else as.character(.data$Description)
        ) |>
        dplyr::select(dplyr::any_of(c("method", result_cols)))

      readr::write_csv(
        result_table,
        file.path(method_dir, sprintf("%s_%s_%s.csv", meta$prefix, file_label, project_name))
      )

    })
  })
}

run_gsea_enrichment <- function(stats, contrasts, genome, out_dir, project_name,
                                methods = c("GO_BP"), parallel_threads = 1,
                                enrichment_lib_dir = NULL,
                                genome_map_file = NULL) {
  methods <- validate_gsea_methods(methods)
  if (length(methods) == 0) {
    return(invisible(NULL))
  }

  purrr::walk(methods, function(method) {
    method_dir <- gsea_csv_dir(out_dir, method)
    unlink(method_dir, recursive = TRUE)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
  })
  message(sprintf("[INFO] GSEA contrasts: %d", length(contrasts)))

  purrr::walk(methods, function(method) {
    meta <- method_output_meta(method, mode = "gsea")
    method_dir <- gsea_csv_dir(out_dir, method)
    method_gene_col <- resolve_enrichment_gene_column(
      genome = genome,
      method = method,
      genome_map_file = genome_map_file
    )
    message(sprintf("[INFO] GSEA %s uses gene column: %s", method, method_gene_col))

    jobs <- lapply(contrasts, function(cv) {
      list(
        method = method,
        numerator = cv[[2]],
        denominator = cv[[3]],
        display_label = format_contrast_label(cv[[2]], cv[[3]])
      )
    })

    run_jobs_parallel(jobs, parallel_threads, function(job) {
      numerator <- job$numerator
      denominator <- job$denominator
      display_label <- job$display_label
      file_label <- format_contrast_file_label(display_label)

      contrast_stats <- stats |>
        dplyr::filter(.data$aspect == display_label)

      gene_list <- build_ranked_gene_list(contrast_stats, gene_col = method_gene_col)
      gsea_table_raw <- run_single_gsea_enrichment(
        method = method,
        gene_list = gene_list,
        genome = genome,
        genome_map_file = genome_map_file,
        enrichment_lib_dir = enrichment_lib_dir
      )

      required_gsea_cols <- c(
        "ID", "Description", "setSize", "enrichmentScore", "NES",
        "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment"
      )
      missing_gsea_cols <- setdiff(required_gsea_cols, colnames(gsea_table_raw))
      if (length(missing_gsea_cols) > 0) {
        for (col_name in missing_gsea_cols) {
          gsea_table_raw[[col_name]] <- NA
        }
      }

      is_kegg_method <- identical(as.character(method), "KEGG")
      gsea_table <- gsea_table_raw |>
        dplyr::mutate(
          method = method,
          contrast = display_label,
          Description = if (is_kegg_method) strip_kegg_species_suffix(.data$Description) else as.character(.data$Description),
          direction = dplyr::case_when(
            !is.na(.data$NES) & .data$NES > 0 ~ numerator,
            !is.na(.data$NES) & .data$NES < 0 ~ denominator,
            TRUE ~ "neutral"
          )
        ) |>
        dplyr::select(dplyr::any_of(c(
          "method", "contrast", "direction", "ID", "Description", "setSize", "enrichmentScore", "NES",
          "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment"
        )))

      readr::write_csv(
        gsea_table,
        file.path(method_dir, sprintf("%s_%s_%s.csv", meta$prefix, file_label, project_name))
      )
    })
  })
}

run_enrichment <- function(stats, genome, out_dir, project_name,
                           methods = c("GO_BP"),
                           enrichment_lib_dir = NULL,
                           genome_map_file = NULL,
                           biomart_lookup_file = NULL,
                           biomart_cache_dir = NULL,
                           biomart_ortholog_dir = NULL,
                           deg_cluster_file = NULL,
                           parallel_threads = 1) {
  stats <- map_stats_to_human_symbols(
    stats = stats,
    genome = genome,
    biomart_lookup_file = biomart_lookup_file,
    biomart_cache_dir = biomart_cache_dir,
    biomart_ortholog_dir = biomart_ortholog_dir
  )
  stats <- stats |>
    dplyr::filter(
      (!is.na(.data$human_gene_symbol) & .data$human_gene_symbol != "") |
        (!is.na(.data$source_gene_symbol) & .data$source_gene_symbol != "")
    )

  contrasts <- build_contrasts_from_stats(stats)
  if (length(contrasts) == 0) {
    stop("No valid contrasts were found in stats$aspect.", call. = FALSE)
  }

  deg_cluster_path <- if (is.null(deg_cluster_file)) "" else as.character(deg_cluster_file)
  if (nzchar(deg_cluster_path) && file.exists(deg_cluster_path)) {
    run_ora_enrichment(
      stats = stats,
      genome = genome,
      out_dir = out_dir,
      project_name = project_name,
      methods = methods,
      parallel_threads = parallel_threads,
      enrichment_lib_dir = enrichment_lib_dir,
      genome_map_file = genome_map_file,
      deg_cluster_file = deg_cluster_file
    )
  } else {
    message(sprintf("[INFO] ORA skipped: deg-clusters file was not found (%s)", deg_cluster_path))
  }

  run_gsea_enrichment(
    stats = stats,
    contrasts = contrasts,
    genome = genome,
    out_dir = out_dir,
    project_name = project_name,
    methods = methods,
    parallel_threads = parallel_threads,
    enrichment_lib_dir = enrichment_lib_dir,
    genome_map_file = genome_map_file
  )
}
