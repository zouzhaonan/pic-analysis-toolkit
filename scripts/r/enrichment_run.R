# 役割:
#   DESeq2 stats から ORA (degPatterns のクラスタ) と GSEA (比較ごと) を回し、CSV を書く。
# 並べ方:
#   method × cluster (ORA) と method × 比較 (GSEA) をすべて 1 つの job の並びにし、
#   重いもの (GSEA、遺伝子セットの多い method) から順に --threads 本ずつ fork で回す。
#   gson は method ごとに親で 1 回だけ読み、子は fork で共有する。
#   fgsea の中の並列 (BiocParallel の既定は全コア) は使わない。job の並列だけで機械を埋める。
# 再現性:
#   GSEA は job ごとに同じ seed (PIC_GSEA_SEED、既定 42) を置いてから回すので、
#   並べ方・--threads・分け方にかかわらず、同じ入力から同じ表が出る。

pic_gsea_seed <- function() {
  s <- suppressWarnings(as.integer(Sys.getenv("PIC_GSEA_SEED", "42")))
  if (is.na(s)) 42L else s
}

# job は list(kind, method, label, cost, run = function() ...)。
# 子の失敗は集めて最後にまとめて止める (ほかの job は最後まで回す)。
# 今のプロセスの常駐メモリ (GB)。Linux は /proc、ほかは ps。読めなければ NA。
pic_self_rss_gb <- function() {
  st <- "/proc/self/status"
  if (file.exists(st)) {
    l <- grep("^VmRSS:", readLines(st, warn = FALSE), value = TRUE)
    if (length(l) == 1) return(as.numeric(sub("^VmRSS:\\s+([0-9]+).*$", "\\1", l)) / 2^20)
  }
  out <- tryCatch(system2("ps", c("-o", "rss=", "-p", Sys.getpid()), stdout = TRUE, stderr = FALSE), error = function(e) "")
  suppressWarnings(as.numeric(trimws(out[[1]])) / 2^20)
}

# worker の数: --threads に従い、減らすのは (1) job の数より多くしない、(2) 親 + worker × 本数 が
# mem_gb を超えない、の 2 つだけ。
# worker 1 本あたりのメモリは、fork したときの親の常駐メモリの 0.5 倍と見積もる。根拠 (pic193_2 mm10、
# 182 job、cgroup の memory.peak): --threads 3 → 3.6 GB、6 → 5.2 GB、8 → 7.1 GB で、1 本増えるごとに
# 約 0.6〜0.7 GB 増え、親 (stats・gson・引き表を読んだあと) は約 1.3 GB。dev035 (hg38) も 8 本で 8.4 GB。
# 子は親の表を書き換えないが、GSEA の途中の作業 (順位表、fgsea、leading edge) と R の GC で私有の
# ページが増えるため、親の大きさに比例させる。
pic_enrich_workers <- function(threads, n_jobs, mem_gb = Inf) {
  n <- max(1L, min(as.integer(threads), as.integer(n_jobs)))
  if (is.finite(mem_gb) && n > 1L) {
    main <- pic_self_rss_gb()
    if (is.na(main)) main <- 1.5
    per <- max(0.25, 0.5 * main)
    fit <- floor((mem_gb - main) / per)
    if (fit < n) {
      message(sprintf("[INFO] enrich: --mem-gb %.1f に収めるため worker を %d から %d に減らします (親 %.1f GB、1 本 %.2f GB と見積もり)",
                      mem_gb, n, max(1L, fit), main, per))
      n <- max(1L, fit)
    }
  }
  n
}

run_jobs_parallel <- function(jobs, parallel_threads, fn = NULL, mem_gb = Inf) {
  if (length(jobs) == 0) {
    return(invisible(NULL))
  }

  threads <- suppressWarnings(as.integer(parallel_threads))
  if (is.na(threads) || threads < 1) {
    threads <- 1L
  }
  threads <- pic_enrich_workers(threads, length(jobs), mem_gb)

  # 子は job の前と後に親が生きているかを見て、親が殺されていたら自分も終わる (analysis_runtime.R)。
  master_pid <- Sys.getpid()
  in_child <- function() !identical(Sys.getpid(), master_pid)
  worker <- function(i) {
    if (in_child()) pic_exit_if_orphaned(master_pid)
    res <- tryCatch(
      {
        if (is.null(fn)) jobs[[i]]$run() else fn(jobs[[i]])
        list(ok = TRUE, idx = i, error = NULL)
      },
      error = function(err) list(ok = FALSE, idx = i, error = conditionMessage(err))
    )
    if (in_child()) pic_exit_if_orphaned(master_pid)
    res
  }

  results <- if (.Platform$OS.type == "unix" && threads > 1 && length(jobs) > 1) {
    # job ごとに fork し、空いた子から次の job を渡す (重さの違う job を偏りなく詰める)。
    parallel::mclapply(
      seq_along(jobs),
      worker,
      mc.cores = min(threads, length(jobs)),
      mc.preschedule = FALSE
    )
  } else {
    lapply(seq_along(jobs), worker)
  }

  failed <- Filter(function(x) !isTRUE(x$ok), results)
  # mclapply は子が落ちると (メモリ不足で殺されたなど) try-error や NULL を返す。
  crashed <- which(vapply(results, function(x) !is.list(x) || is.null(x$ok), logical(1)))
  for (i in crashed) {
    failed[[length(failed) + 1]] <- list(ok = FALSE, idx = i, error = paste(as.character(results[[i]]), collapse = " "))
  }
  failed <- Filter(function(x) is.list(x) && !isTRUE(x$ok), failed)
  if (length(failed) > 0) {
    details <- vapply(
      failed,
      function(x) {
        job <- jobs[[x$idx]]
        sprintf(
          "%s method=%s label=%s error=%s",
          as.character(job$kind %||% ""),
          as.character(job$method),
          as.character(job$label %||% sprintf("job_%d", as.integer(x$idx))),
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

`%||%` <- function(a, b) if (is.null(a)) b else a

# 書きかけの CSV を残さない (止められたときに、レポートが半端な表を読まないように)。
write_csv_atomic <- function(tbl, path) {
  tmp <- paste0(path, ".tmp", Sys.getpid())
  readr::write_csv(tbl, tmp)
  file.rename(tmp, path)
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

# method_dir を作り、この project の旧 CSV (*_<project_name>.csv) だけを消す。
# pic all は全 genome が同じ out_dir/enrich を共有するので、フォルダごと消すと
# 先に走った genome の結果まで消えてしまう。
reset_method_dir <- function(method_dir, project_name) {
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
  old <- list.files(method_dir, full.names = TRUE)
  old <- old[endsWith(basename(old), sprintf("_%s.csv", project_name))]
  unlink(old)
}


# geneID は遺伝子の symbol を "/" で並べたもの。そのすぐ右に、同じ並びの ens_gene (geneID_ens_gene) を置く
# (gene_ids: symbol -> ens_gene。load_precomputed_deg_clusters の gene_ids)。
ora_result_table <- function(method, cluster_id, cluster_genes, universe_genes, gson_obj, anno_env = NULL, gene_ids = NULL) {
  display_label <- cluster_id
  if (length(cluster_genes) == 0) {
    empty_cols <- c(
      "method", "cluster", "ID", "Description", "GeneRatio", "BgRatio",
      "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "geneID_ens_gene", "Count"
    )
    empty_tbl <- as.data.frame(setNames(rep(list(character()), length(empty_cols)), empty_cols), stringsAsFactors = FALSE)
    return(tibble::as_tibble(empty_tbl))
  }

  cluster_res <- safe_as_tibble(run_local_ora(cluster_genes, universe_genes, gson_obj, anno_env)) |>
    dplyr::mutate(direction = cluster_id)

  required_ora_cols <- c(
    "ID", "Description", "GeneRatio", "BgRatio",
    "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "direction"
  )
  numeric_ora_cols <- c("RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "Count")
  missing <- setdiff(required_ora_cols, colnames(cluster_res))
  for (col_name in missing) {
    if (col_name %in% numeric_ora_cols) {
      cluster_res[[col_name]] <- NA_real_
    } else {
      cluster_res[[col_name]] <- NA_character_
    }
  }
  for (col_name in required_ora_cols) {
    if (col_name %in% numeric_ora_cols) {
      cluster_res[[col_name]] <- suppressWarnings(as.numeric(cluster_res[[col_name]]))
    } else {
      cluster_res[[col_name]] <- as.character(cluster_res[[col_name]])
    }
  }

  result_cols <- c(
    "cluster", "ID", "Description", "GeneRatio", "BgRatio",
    "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "geneID_ens_gene", "Count"
  )
  is_kegg_method <- identical(as.character(method), "KEGG")
  cluster_res |>
    dplyr::filter(!is.na(.data$pvalue)) |>
    dplyr::arrange(.data$pvalue) |>
    dplyr::mutate(
      cluster = display_label,
      method = method,
      Description = if (is_kegg_method) strip_kegg_species_suffix(.data$Description) else as.character(.data$Description),
      geneID_ens_gene = pic_genes_to_ids(.data$geneID, gene_ids)
    ) |>
    dplyr::select(dplyr::any_of(c("method", result_cols)))
}

# core_enrichment は遺伝子の symbol を "/" で並べたもの。そのすぐ右に、同じ並びの ens_gene
# (core_enrichment_ens_gene) を置く (gene_ids: symbol -> 順位に使った行の ens_gene。ranked_gene_table)。
gsea_result_table <- function(method, numerator, denominator, display_label, gene_list, gson_obj, gene_ids = NULL) {
  gsea_table_raw <- safe_as_tibble(run_local_gsea(gene_list, gson_obj, seed = pic_gsea_seed()))

  required_gsea_cols <- c(
    "ID", "Description", "setSize", "enrichmentScore", "NES",
    "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment"
  )
  for (col_name in setdiff(required_gsea_cols, colnames(gsea_table_raw))) {
    gsea_table_raw[[col_name]] <- NA
  }

  is_kegg_method <- identical(as.character(method), "KEGG")
  gsea_table_raw |>
    dplyr::mutate(
      method = method,
      contrast = display_label,
      Description = if (is_kegg_method) strip_kegg_species_suffix(.data$Description) else as.character(.data$Description),
      direction = dplyr::case_when(
        !is.na(.data$NES) & .data$NES > 0 ~ numerator,
        !is.na(.data$NES) & .data$NES < 0 ~ denominator,
        TRUE ~ "neutral"
      ),
      core_enrichment_ens_gene = pic_genes_to_ids(.data$core_enrichment, gene_ids)
    ) |>
    dplyr::select(dplyr::any_of(c(
      "method", "contrast", "direction", "ID", "Description", "setSize", "enrichmentScore", "NES",
      "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment", "core_enrichment_ens_gene"
    )))
}

# gson の中で minGSSize..maxGSSize に入る遺伝子セットの数 (job の重さの目安)。
gson_set_count <- function(gson_obj) {
  if (is.null(gson_obj)) return(0L)
  g2 <- tryCatch(gson_obj@gsid2gene, error = function(e) NULL)
  if (is.null(g2) || nrow(g2) == 0) return(0L)
  n <- table(g2[[1]])
  sum(n >= 10 & n <= 500)
}

# ORA の job (method × cluster) を作る。method の gson と universe は親で 1 回だけ用意する。
build_ora_jobs <- function(stats, genome, out_dir, project_name, methods,
                           genome_map_file, deg_cluster_file, gson_cache) {
  methods <- validate_ora_methods(methods)
  # ORA の引き表は method ごとに親で 1 回だけ作る (子は fork で共有する)。
  anno_cache <- stats::setNames(lapply(methods, function(m) gson_to_anno_env(gson_cache[[m]])), methods)
  if (length(methods) == 0) return(list())
  if (is.null(deg_cluster_file) || !nzchar(as.character(deg_cluster_file))) {
    stop("ORA requires --deg-clusters.", call. = FALSE)
  }

  for (method in methods) reset_method_dir(ora_csv_dir(out_dir, method), project_name)

  first_clusters <- load_precomputed_deg_clusters(
    cluster_file = deg_cluster_file,
    stats = stats,
    gene_col = resolve_enrichment_gene_column(genome = genome, method = methods[[1]], genome_map_file = genome_map_file)
  )
  message(sprintf("[INFO] ORA clustering by degPatterns: %d clusters", nrow(first_clusters)))

  jobs <- list()
  for (method in methods) {
    method_gene_col <- resolve_enrichment_gene_column(genome = genome, method = method, genome_map_file = genome_map_file)
    message(sprintf("[INFO] ORA %s uses gene column: %s", method, method_gene_col))

    clusters <- load_precomputed_deg_clusters(cluster_file = deg_cluster_file, stats = stats, gene_col = method_gene_col)
    if (nrow(clusters) == 0) {
      message(sprintf("[INFO] ORA %s skipped: no clusters after mapping", method))
      next
    }
    universe_genes <- stats |>
      dplyr::mutate(enrich_gene_symbol = as.character(.data[[method_gene_col]])) |>
      dplyr::filter(!is.na(.data$enrich_gene_symbol), .data$enrich_gene_symbol != "") |>
      dplyr::pull(.data$enrich_gene_symbol) |>
      unique()
    meta <- method_output_meta(method, mode = "ora")
    method_dir <- ora_csv_dir(out_dir, method)

    for (idx in seq_len(nrow(clusters))) {
      # local の中で値を束ねる (ループの変数は回るたびに上書きされるので、job から直に参照しない)。
      local({
        m <- method
        univ <- universe_genes
        cluster_id <- as.character(clusters$cluster_id[[idx]])
        cluster_genes <- unique(as.character(clusters$genes[[idx]]))
        cluster_ids <- clusters$gene_ids[[idx]]
        out_csv <- file.path(method_dir, sprintf("%s_%s_%s.csv", meta$prefix, format_contrast_file_label(cluster_id), project_name))
        jobs[[length(jobs) + 1]] <<- list(
          kind = "ORA", method = m, label = cluster_id, cost = 1,
          run = function() {
            write_csv_atomic(ora_result_table(m, cluster_id, cluster_genes, univ, gson_cache[[m]], anno_cache[[m]],
                                              gene_ids = cluster_ids), out_csv)
          }
        )
      })
    }
  }
  jobs
}

# GSEA の job (method × 比較) を作る。
build_gsea_jobs <- function(stats, contrasts, genome, out_dir, project_name, methods,
                            genome_map_file, gson_cache) {
  methods <- validate_gsea_methods(methods)
  if (length(methods) == 0) return(list())

  for (method in methods) reset_method_dir(gsea_csv_dir(out_dir, method), project_name)
  message(sprintf("[INFO] GSEA contrasts: %d", length(contrasts)))

  jobs <- list()
  for (method in methods) {
    method_gene_col <- resolve_enrichment_gene_column(genome = genome, method = method, genome_map_file = genome_map_file)
    message(sprintf("[INFO] GSEA %s uses gene column: %s", method, method_gene_col))
    meta <- method_output_meta(method, mode = "gsea")
    method_dir <- gsea_csv_dir(out_dir, method)
    cost <- max(1, gson_set_count(gson_cache[[method]]))

    for (cv in contrasts) {
      local({
        m <- method
        gene_col <- method_gene_col
        numerator <- cv[[2]]
        denominator <- cv[[3]]
        display_label <- format_contrast_label(numerator, denominator)
        out_csv <- file.path(method_dir, sprintf("%s_%s_%s.csv", meta$prefix, format_contrast_file_label(display_label), project_name))
        jobs[[length(jobs) + 1]] <<- list(
          kind = "GSEA", method = m, label = display_label, cost = cost,
          run = function() {
            contrast_stats <- stats |>
              dplyr::filter(.data$aspect == display_label)
            rank_tbl <- ranked_gene_table(contrast_stats, gene_col = gene_col)
            gene_list <- ranked_list_from_table(rank_tbl)
            gene_ids <- stats::setNames(as.character(rank_tbl$ens_gene), rank_tbl$gene_id)
            write_csv_atomic(gsea_result_table(m, numerator, denominator, display_label, gene_list, gson_cache[[m]],
                                               gene_ids = gene_ids), out_csv)
          }
        )
      })
    }
  }
  jobs
}

run_enrichment <- function(stats, genome, out_dir, project_name,
                           methods = c("GO_BP"),
                           enrichment_lib_dir = NULL,
                           genome_map_file = NULL,
                           biomart_lookup_file = NULL,
                           biomart_cache_dir = NULL,
                           biomart_ortholog_dir = NULL,
                           deg_cluster_file = NULL,
                           parallel_threads = 1,
                           analyses = c("ORA", "GSEA"),
                           mem_gb = Inf) {
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

  # GSEA の leading edge を同じ値の速い版にする (DOSE 4.0.0 のときだけ。enrichment_gsea_fast.R)。
  if ("GSEA" %in% analyses && !pic_gsea_fast_install()) {
    message("[INFO] GSEA leading edge: DOSE の元の計算を使います")
  }

  # gson は method ごとに 1 回だけ読む (子は fork で共有する)。
  all_methods <- unique(c(validate_ora_methods(methods), validate_gsea_methods(methods)))
  gson_cache <- stats::setNames(lapply(all_methods, function(m) {
    read_local_gson(genome = genome, method = m, genome_map_file = genome_map_file, enrichment_lib_dir = enrichment_lib_dir)
  }), all_methods)

  jobs <- list()
  deg_cluster_path <- if (is.null(deg_cluster_file)) "" else as.character(deg_cluster_file)
  if ("ORA" %in% analyses) {
    if (nzchar(deg_cluster_path) && file.exists(deg_cluster_path)) {
      jobs <- c(jobs, build_ora_jobs(
        stats = stats, genome = genome, out_dir = out_dir, project_name = project_name,
        methods = methods, genome_map_file = genome_map_file,
        deg_cluster_file = deg_cluster_file, gson_cache = gson_cache
      ))
    } else {
      message(sprintf("[INFO] ORA skipped: deg-clusters file was not found (%s)", deg_cluster_path))
    }
  }
  if ("GSEA" %in% analyses) {
    jobs <- c(jobs, build_gsea_jobs(
      stats = stats, contrasts = contrasts, genome = genome, out_dir = out_dir,
      project_name = project_name, methods = methods,
      genome_map_file = genome_map_file, gson_cache = gson_cache
    ))
  }

  # 重いものから回す (最後に重い job が 1 本だけ残るのを避ける)。GSEA は ORA よりずっと重い。
  if (length(jobs) > 0) {
    weight <- vapply(jobs, function(j) if (identical(j$kind, "GSEA")) 1e6 + j$cost else j$cost, numeric(1))
    jobs <- jobs[order(-weight, seq_along(jobs))]
  }
  n_gsea <- sum(vapply(jobs, function(j) identical(j$kind, "GSEA"), logical(1)))
  message(sprintf("[INFO] Enrichment jobs: %d (GSEA %d, ORA %d), threads=%d",
                  length(jobs), n_gsea, length(jobs) - n_gsea, as.integer(parallel_threads)))
  run_jobs_parallel(jobs, parallel_threads, mem_gb = mem_gb)
}
