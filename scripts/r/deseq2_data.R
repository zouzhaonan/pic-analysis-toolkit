load_umi_tables <- function(def, count_dir) {
  if (!dir.exists(count_dir)) {
    stop(sprintf("Count directory was not found: %s", count_dir), call. = FALSE)
  }

  normalize_count_table <- function(tbl, path_label) {
    nm <- tolower(colnames(tbl))
    gene_idx <- which(nm %in% c("ens_gene", "gene", "feature", "gene_id"))
    cell_idx <- which(nm %in% c("barcode", "cell", "cell_barcode"))
    count_idx <- which(nm %in% c("count", "counts", "umi", "umi_count"))

    if (length(gene_idx) == 0 || length(cell_idx) == 0 || length(count_idx) == 0) {
      stop(
        sprintf(
          "Invalid count table format: %s (required columns for gene/cell/count were not found)",
          path_label
        ),
        call. = FALSE
      )
    }

    tibble::tibble(
      ens_gene = as.character(tbl[[gene_idx[[1]]]]),
      barcode = as.character(tbl[[cell_idx[[1]]]]),
      count = suppressWarnings(as.numeric(tbl[[count_idx[[1]]]]))
    ) |>
      dplyr::filter(!is.na(.data$ens_gene), .data$ens_gene != "", !is.na(.data$barcode), .data$barcode != "")
  }

  umi <- def |>
    dplyr::distinct(.data$count_prefix) |>
    dplyr::mutate(
      count_file = file.path(count_dir, paste0(.data$count_prefix, ".txt.gz")),
      data = purrr::map(
        .data$count_file,
        ~ {
          raw_tbl <- suppressMessages(readr::read_tsv(
            .x,
            show_col_types = FALSE,
            progress = FALSE
          ))
          normalize_count_table(raw_tbl, .x)
        }
      )
    ) |>
    tidyr::unnest("data")

  umi |>
    dplyr::inner_join(
      dplyr::select(def, "count_prefix", "barcode", "sample"),
      by = c("count_prefix", "barcode")
    ) |>
    dplyr::select(-"count_prefix", -"barcode")
}

build_count_matrix <- function(umi, e2g) {
  # 遺伝子 × サンプルに足し上げて横に広げる。data.table で行う (dplyr の group_by + pivot_wider と
  # 同じ表: 行は ens_gene を C ロケールで並べた順、列はその並びで最初に出てくるサンプルの順、無い組は 0)。
  dt <- data.table::data.table(
    ens_gene = as.character(umi$ens_gene),
    sample = as.character(umi$sample),
    count = as.numeric(umi$count)
  )
  dt <- dt[!is.na(ens_gene) & ens_gene != "" & !is.na(sample) & sample != ""]
  agg <- dt[, list(count = sum(count, na.rm = TRUE)), keyby = c("ens_gene", "sample")]
  if (nrow(agg) == 0) {
    sample_mat <- tibble::tibble(ens_gene = character())
  } else {
    sample_cols <- unique(agg$sample)
    wide <- data.table::dcast(agg, ens_gene ~ sample, value.var = "count", fill = 0)
    data.table::setcolorder(wide, c("ens_gene", sample_cols))
    sample_mat <- tibble::as_tibble(as.data.frame(wide, stringsAsFactors = FALSE))
  }

  if (nrow(sample_mat) == 0) {
    return(tibble::tibble(
      ens_gene = character(),
      ext_gene = character(),
      biotype = character(),
      chr = character()
    ))
  }

  if (nrow(e2g) == 0) {
    return(
      sample_mat |>
        dplyr::mutate(
          ext_gene = .data$ens_gene,
          biotype = NA_character_,
          chr = NA_character_,
          .after = "ens_gene"
        )
    )
  }

  # ens_gene ごとに、空でない最初の ext_gene / biotype / chr (左結合に使うだけなので行の順は問わない)。
  e <- data.table::data.table(
    ens_gene = as.character(e2g$ens_gene),
    ext_gene = as.character(e2g$ext_gene),
    biotype = as.character(e2g$biotype),
    chr = as.character(e2g$chr)
  )
  e <- e[!is.na(ens_gene) & ens_gene != ""]
  ids <- unique(e$ens_gene)
  first_nonempty <- function(col) {
    v <- e[[col]]
    ok <- !is.na(v) & v != ""
    hit <- e$ens_gene[ok]
    keep <- !duplicated(hit)
    v[ok][keep][match(ids, hit[keep])]
  }
  e2g_uniq <- tibble::tibble(
    ens_gene = ids,
    ext_gene = first_nonempty("ext_gene"),
    biotype = first_nonempty("biotype"),
    chr = first_nonempty("chr")
  )

  sample_mat |>
    dplyr::left_join(e2g_uniq, by = "ens_gene") |>
    dplyr::mutate(ext_gene = dplyr::coalesce(.data$ext_gene, .data$ens_gene)) |>
    dplyr::relocate(
      dplyr::all_of(c("ext_gene", "biotype", "chr")),
      .after = "ens_gene"
    )
}

format_group_label <- function(group_name) {
  label <- as.character(group_name)
  label <- tolower(label)
  gsub("[[:space:]]+", "_", label)
}

format_contrast_label <- function(numerator, denominator) {
  paste0(format_group_label(numerator), " / ", format_group_label(denominator))
}

# サイズ因子を推定して dds に設定する。
#
# 既定は poscounts (ゼロが多いカウントに対する DESeq2 の推奨)。ただし UMI カウントが
# 極端に疎な場合 (大半の遺伝子が 1-2 カウント) に破綻する。深いライブラリは
# 「遺伝子あたりのカウントが増える」のではなく「検出遺伝子数が増える」形で深くなるため、
# 「正の値を持つ遺伝子における比の中央値」が全サンプルで 1 になり、深度差を検出できない。
# その状態では normalizedCountTable が生カウントと同一になり、深いサンプルほど
# 高発現に見えるという系統的なバイアスが全遺伝子に入る。
# ここでは退化を検知し、ライブラリサイズ正規化にフォールバックする。
pic_set_size_factors <- function(dds) {
  lib <- colSums(DESeq2::counts(dds))
  lib <- pmax(lib, 1)

  sf <- tryCatch(
    DESeq2::sizeFactors(DESeq2::estimateSizeFactors(dds, type = "poscounts")),
    error = function(err) NULL
  )

  degenerate <- is.null(sf) || any(!is.finite(sf)) || any(sf <= 0)
  if (!degenerate && length(sf) >= 2) {
    # サイズ因子が深度差をどれだけ追随できているかで判定する。
    # 正常なら log(sizeFactor) は log(ライブラリサイズ) と同程度にばらつく。
    # 退化していると、深度が何倍違ってもサイズ因子はほぼ動かない。
    # (異常サンプル 1 本だけが外れる場合に単純な最大/最小比では検知できないため、
    #  ばらつきの比で判定する。)
    sd_sf <- stats::sd(log(sf))
    sd_lib <- stats::sd(log(lib))
    degenerate <- is.finite(sd_lib) && sd_lib > 0 &&
      (sd_sf < 1e-8 || (sd_lib > 0.1 && sd_sf < 0.2 * sd_lib))
  }

  if (degenerate) {
    sf_ratio <- if (is.null(sf)) NA_real_ else max(sf) / min(sf)
    message(
      "[WARN] poscounts のサイズ因子がライブラリ深度を反映していません ",
      sprintf("(sizeFactor 比 %.2f / ライブラリサイズ比 %.2f)。",
              sf_ratio, max(lib) / min(lib)),
      "ライブラリサイズ正規化にフォールバックします。"
    )
    sf <- lib / exp(mean(log(lib)))
  }

  # 実際に採用した方法を記録する (HTML レポートの Materials & Methods が参照する)。
  assign("pic_size_factor_method",
         if (degenerate) "libsize" else "poscounts",
         envir = globalenv())

  DESeq2::sizeFactors(dds) <- sf
  dds
}

# summary/analysis_params.tsv に key<TAB>value を追記/更新する。
pic_write_analysis_param <- function(out_dir, key, value) {
  f <- file.path(dirname(normalizePath(out_dir, mustWork = FALSE)),
                 "summary", "analysis_params.tsv")
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  rows <- if (file.exists(f)) {
    tryCatch(utils::read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
  } else NULL
  if (is.null(rows) || !all(c("key", "value") %in% names(rows))) {
    rows <- data.frame(key = character(0), value = character(0), stringsAsFactors = FALSE)
  }
  rows <- rows[rows$key != key, , drop = FALSE]
  rows <- rbind(rows, data.frame(key = key, value = as.character(value),
                                 stringsAsFactors = FALSE))
  tmp <- pic_tmp_path(f)   # 一時ファイル -> rename (書きかけを残さない。deseq2_common.R)
  utils::write.table(rows, tmp, sep = "\t", quote = FALSE, row.names = FALSE)
  pic_commit_tmp(tmp, f)
}

# カウント行列 (UMI_count の表) から DESeqDataSet を作り、サイズ因子を確定させる。
# fit と rlog の両方がこれから始める (同じ入力なら同じ dds・同じサイズ因子になる)。
build_deseq_dds <- function(mat, def) {
  mat_filtered <- mat |>
    dplyr::filter(.data$chr != "MT" | is.na(.data$chr)) |>
    dplyr::filter(dplyr::if_any(-(1:4), ~ .x > 0))

  count_matrix <- as.matrix(dplyr::select(mat_filtered, -(1:4)))
  if (ncol(count_matrix) == 0 || nrow(count_matrix) == 0) {
    stop(
      "Count matrix is empty after merging deftable and count files. Check count_prefix/barcode consistency.",
      call. = FALSE
    )
  }
  mode(count_matrix) <- "numeric"
  rownames(count_matrix) <- mat_filtered$ens_gene

  label <- def |>
    dplyr::distinct(.data$sample, .data$group) |>
    dplyr::filter(.data$sample %in% colnames(count_matrix))

  label$group <- factor(label$group, levels = unique(label$group))
  label_df <- as.data.frame(label)
  rownames(label_df) <- label_df$sample

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix[, label_df$sample, drop = FALSE],
    colData = label_df,
    design = ~group
  )
  # サイズ因子を先に確定させる (退化時はライブラリサイズにフォールバック)。
  # 設定済みなら DESeq() は再推定しないため sfType は渡さない。
  dds <- pic_set_size_factors(dds)

  list(dds = dds, label = label_df)
}

# rlog (blind = TRUE) を fork した子で回し始める (--part all で threads > 1 のとき、DESeq() と同時に)。
# rlog は design を ~1 にして分散を推定し直すので、DESeq() の前の dds から作っても後の dds から
# 作っても同じ行列になる。返り値を collect_rlog_matrix() に渡して受け取る。
start_rlog_job <- function(dds, master_pid = Sys.getpid()) {
  force(master_pid)   # 既定値の評価を fork の前に (子の中で評価すると子の pid になる)
  parallel::mcparallel({
    pic_start_child_guard(master_pid)   # 親が殺されたら、rlog の途中でもこの子を止める (deseq2_common.R)
    r <- compute_rlog_matrix(dds)
    pic_exit_if_orphaned(master_pid)    # 親が殺されていたら結果を送らずに終わる (analysis_runtime.R)
    r
  }, silent = TRUE)
}

# DESeq() と比較ごとの results()。
#   - results() は比較ごとに fork して並べる (比較どうしは独立。1 比較 1〜2 秒かかる)。
#     worker の数は pic_n_workers (threads・比較の数・メモリ) で決める。busy は同時に走っている
#     ほかの子 (rlog) の分: threads と mem_gb からその分を引く。
#   - DESeq() 自体の parallel = TRUE は使わない (36 サンプルで試すと、子に dds を渡す手間で 1 本より遅い)。
fit_deseq <- function(dds, contrasts, threads = 1L, mem_gb = Inf, master_pid = Sys.getpid(),
                      busy_threads = 0L, busy_gb = 0) {
  force(master_pid)   # 既定値の評価を fork の前に (子の中で評価すると子の pid になる)
  deseq_fit <- tryCatch(
    DESeq2::DESeq(dds),
    error = function(err) err
  )

  if (inherits(deseq_fit, "error")) {
    no_rep_msg <- grepl(
      "estimation of dispersion is not possible",
      conditionMessage(deseq_fit),
      fixed = TRUE
    )
    if (!isTRUE(no_rep_msg)) {
      stop(deseq_fit)
    }

    nc <- DESeq2::counts(dds, normalized = TRUE)
    rng <- apply(log2(nc + 1), 1, function(x) diff(range(x)))
    ok <- rowMeans(nc) > 5
    stable <- ok & rng <= stats::quantile(rng[ok], 0.5)
    ref <- dds[stable, ]
    DESeq2::design(ref) <- ~1
    ref <- DESeq2::estimateDispersionsFit(DESeq2::estimateDispersionsGeneEst(ref))
    DESeq2::dispersions(dds) <- pmax(DESeq2::dispersionFunction(ref)(rowMeans(nc)), 1e-8)
    dds <- DESeq2::nbinomWaldTest(dds)
    assign("pic_dispersion_method", "stable_trend", envir = globalenv())
  } else {
    dds <- deseq_fit
    assign("pic_dispersion_method", "group", envir = globalenv())
  }

  one_result <- function(cv) {
    DESeq2::results(
      dds,
      contrast = cv,
      lfcThreshold = log2(1),
      independentFiltering = FALSE,
      cooksCutoff = FALSE
    )
  }
  mem <- deseq2_mem_estimate_gb(nrow(dds), ncol(dds))
  workers <- pic_n_workers(max(1L, threads - busy_threads), length(contrasts),
                           per_worker_gb = mem$results_worker, mem_gb = mem_gb - busy_gb, main_gb = mem$main)
  results_list <- if (workers > 1L && .Platform$OS.type == "unix") {
    message(sprintf("[INFO] results(): %d contrasts on %d workers", length(contrasts), workers))
    res <- parallel::mclapply(contrasts, function(cv) {
      pic_exit_if_orphaned(master_pid)
      r <- one_result(cv)
      pic_exit_if_orphaned(master_pid)
      r
    }, mc.cores = workers)
    bad <- vapply(res, function(r) !inherits(r, "DESeqResults"), logical(1))
    if (any(bad)) stop(sprintf("DESeq2::results failed for: %s (%s)", paste(names(contrasts)[bad], collapse = ", "),
                               paste(unique(vapply(res[bad], function(r) paste(as.character(r), collapse = " "), character(1))), collapse = " | ")),
                       call. = FALSE)
    res
  } else {
    lapply(contrasts, one_result)
  }

  list(dds = dds, stats = results_list)
}

build_stats_tables <- function(results_list, contrasts, e2g, fdr) {
  annotate <- function(df) {
    if (nrow(e2g) == 0) {
      return(df)
    }
    dplyr::left_join(
      df,
      dplyr::distinct(
        e2g,
        .data$ens_gene,
        .data$ext_gene,
        .data$biotype,
        .data$chr
      ),
      by = "ens_gene"
    )
  }

  stats <- purrr::imap(results_list, function(res, contrast_name) {
    contrast_values <- contrasts[[contrast_name]]
    display_label <- format_contrast_label(
      contrast_values[[2]],
      contrast_values[[3]]
    )

    tibble::as_tibble(res, rownames = "ens_gene") |>
      dplyr::mutate(aspect = display_label, .before = 1) |>
      annotate() |>
      dplyr::relocate(
        dplyr::any_of(c("ens_gene", "ext_gene", "biotype", "chr")),
        .after = "aspect"
      )
  }) |>
    dplyr::bind_rows() |>
    dplyr::arrange(is.na(.data$pvalue), .data$pvalue)

  list(stats = stats)
}

build_normalized_counts <- function(dds, e2g) {
  normcount <- tibble::as_tibble(
    DESeq2::counts(dds, normalized = TRUE),
    rownames = "ens_gene"
  )

  if (nrow(e2g) > 0) {
    normcount <- dplyr::left_join(
      normcount,
      dplyr::distinct(
        e2g,
        .data$ens_gene,
        .data$ext_gene,
        .data$biotype,
        .data$chr
      ),
      by = "ens_gene"
    )
    normcount <- dplyr::relocate(
      normcount,
      dplyr::any_of(c("ext_gene", "biotype", "chr")),
      .after = "ens_gene"
    )
  }

  normcount
}

build_deg_count_summary <- function(stats, contrasts, fdr) {
  contrast_names <- names(contrasts)
  if (length(contrast_names) == 0) {
    return(tibble::tibble())
  }

  purrr::imap_dfr(contrasts, function(contrast_values, contrast_name) {
    display_label <- format_contrast_label(
      contrast_values[[2]],
      contrast_values[[3]]
    )
    up_count <- sum(
      stats$aspect == display_label &
        stats$padj < fdr &
        stats$log2FoldChange > 0,
      na.rm = TRUE
    )
    down_count <- sum(
      stats$aspect == display_label &
        stats$padj < fdr &
        stats$log2FoldChange < 0,
      na.rm = TRUE
    )

    tibble::tibble(
      contrast = c(display_label, display_label),
      biased_group = c(contrast_values[[2]], contrast_values[[3]]),
      deg_count = c(up_count, down_count)
    )
  })
}

# rlog (blind = TRUE) の行列。DEG のクラスタ分けと PCA が同じものを使うので、main で 1 回だけ
# 計算して両方に渡す (rlog はサンプルが多いと DESeq2 の段で一番重い)。失敗したら NULL。
compute_rlog_matrix <- function(dds) {
  tryCatch(
    SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE)),
    error = function(err) NULL
  )
}

# start_rlog_job が子で回した rlog を受け取る。子が無い・落ちたときはここで計算する (同じ行列)。
collect_rlog_matrix <- function(dds, rlog_job = NULL) {
  if (!is.null(rlog_job)) {
    res <- parallel::mccollect(rlog_job, wait = TRUE)[[1]]
    if (is.null(res) || is.matrix(res)) return(res)
    message("[WARN] rlog の子が失敗したので、計算し直します")
  }
  compute_rlog_matrix(dds)
}

build_rlog_matrix_for_degpatterns <- function(dds, rlog_mat = NULL) {
  rlog_try <- if (is.null(rlog_mat)) compute_rlog_matrix(dds) else rlog_mat
  if (!is.null(rlog_try)) {
    return(rlog_try)
  }

  norm <- DESeq2::counts(dds, normalized = TRUE)
  log2(norm + 1)
}

build_degpattern_outputs <- function(
  dds,
  label,
  stats,
  e2g,
  fdr,
  min_cluster_size = 5L,
  rlog_mat = NULL
) {
  deg_genes <- stats |>
    dplyr::filter(
      !is.na(.data$ens_gene),
      .data$ens_gene != "",
      !is.na(.data$padj),
      .data$padj < fdr
    ) |>
    dplyr::pull(.data$ens_gene) |>
    unique()

  if (length(deg_genes) == 0) {
    return(list(
      cluster_gene = tibble::tibble(),
      cluster_summary = tibble::tibble(),
      cluster_profile = tibble::tibble(),
      plot = NULL
    ))
  }

  rlog_mat <- build_rlog_matrix_for_degpatterns(dds, rlog_mat)
  deg_genes <- intersect(deg_genes, rownames(rlog_mat))
  if (length(deg_genes) == 0) {
    return(list(
      cluster_gene = tibble::tibble(),
      cluster_summary = tibble::tibble(),
      cluster_profile = tibble::tibble(),
      plot = NULL
    ))
  }

  mat <- rlog_mat[deg_genes, , drop = FALSE]
  sample_meta <- label |>
    dplyr::select("sample", "group") |>
    dplyr::distinct()
  sample_meta <- as.data.frame(sample_meta, stringsAsFactors = FALSE)
  rownames(sample_meta) <- sample_meta$sample
  available_samples <- intersect(colnames(mat), rownames(sample_meta))
  mat <- mat[, available_samples, drop = FALSE]
  sample_meta <- sample_meta[available_samples, , drop = FALSE]
  sample_meta$group <- as.character(sample_meta$group)
  if (ncol(mat) < 2 || nrow(sample_meta) < 2) {
    return(list(
      cluster_gene = tibble::tibble(),
      cluster_summary = tibble::tibble(),
      cluster_profile = tibble::tibble(),
      plot = NULL
    ))
  }

  group_names <- unique(sample_meta$group)
  group_profile <- vapply(
    group_names,
    function(grp) {
      cols <- rownames(sample_meta)[sample_meta$group == grp]
      rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
    },
    numeric(nrow(mat))
  )
  # vapply は FUN.VALUE の長さが 1 (= DEG がちょうど 1 遺伝子) のとき matrix ではなく
  # vector を返すため、そのままでは次元が 遺伝子 × 群 にならない。明示的に整形する。
  group_profile <- matrix(
    group_profile,
    nrow = nrow(mat),
    ncol = length(group_names),
    dimnames = list(rownames(mat), group_names)
  )

  if (nrow(group_profile) == 0 || ncol(group_profile) == 0) {
    return(list(
      cluster_gene = tibble::tibble(),
      cluster_summary = tibble::tibble(),
      cluster_profile = tibble::tibble(),
      plot = NULL
    ))
  }

  cluster_ids <- rep(1L, nrow(group_profile))
  names(cluster_ids) <- rownames(group_profile)
  if (nrow(group_profile) >= 2 && ncol(group_profile) >= 2) {
    cor_mat <- suppressWarnings(stats::cor(
      t(group_profile),
      use = "pairwise.complete.obs",
      method = "pearson"
    ))
    cor_mat[is.na(cor_mat)] <- 0
    diag(cor_mat) <- 1
    hc <- stats::hclust(stats::as.dist(1 - cor_mat), method = "average")
    max_k <- max(1L, floor(nrow(group_profile) / max(1L, as.integer(min_cluster_size))))
    k <- min(max_k, 12L)
    k <- max(k, 1L)
    cluster_ids <- stats::cutree(hc, k = k)
  }

  cluster_tbl <- tibble::tibble(
    ens_gene = names(cluster_ids),
    cluster_num = as.integer(cluster_ids)
  )
  cluster_size_tbl <- cluster_tbl |>
    dplyr::count(.data$cluster_num, name = "n_genes")
  valid_clusters <- cluster_size_tbl |>
    dplyr::filter(.data$n_genes >= as.integer(min_cluster_size)) |>
    dplyr::pull(.data$cluster_num)
  if (length(valid_clusters) > 0) {
    cluster_tbl <- cluster_tbl |>
      dplyr::filter(.data$cluster_num %in% valid_clusters)
  } else {
    cluster_tbl <- cluster_tbl
  }

  if (nrow(cluster_tbl) == 0) {
    return(list(
      cluster_gene = tibble::tibble(),
      cluster_summary = tibble::tibble(),
      cluster_profile = tibble::tibble(),
      plot = NULL
    ))
  }

  cluster_gene <- cluster_tbl |>
    dplyr::transmute(
      ens_gene = as.character(.data$ens_gene),
      cluster_id = paste0("cluster_", as.character(.data$cluster_num))
    ) |>
    dplyr::filter(
      !is.na(.data$ens_gene),
      .data$ens_gene != "",
      !is.na(.data$cluster_id),
      .data$cluster_id != ""
    ) |>
    dplyr::distinct(.data$cluster_id, .data$ens_gene)

  if (nrow(cluster_gene) == 0) {
    return(list(
      cluster_gene = tibble::tibble(),
      cluster_summary = tibble::tibble(),
      cluster_profile = tibble::tibble(),
      plot = NULL
    ))
  }

  if (nrow(e2g) > 0) {
    cluster_gene <- cluster_gene |>
      dplyr::left_join(
        dplyr::distinct(
          e2g,
          .data$ens_gene,
          .data$ext_gene,
          .data$biotype,
          .data$chr
        ),
        by = "ens_gene"
      ) |>
      dplyr::relocate(
        dplyr::any_of(c("ext_gene", "biotype", "chr")),
        .after = "ens_gene"
      )
  }

  profile_long <- tibble::as_tibble(mat, rownames = "ens_gene") |>
    tidyr::pivot_longer(
      cols = -"ens_gene",
      names_to = "sample",
      values_to = "rlog_expr"
    ) |>
    dplyr::left_join(
      cluster_gene |> dplyr::select("cluster_id", "ens_gene"),
      by = "ens_gene"
    ) |>
    dplyr::left_join(tibble::as_tibble(sample_meta), by = "sample") |>
    dplyr::filter(!is.na(.data$cluster_id), !is.na(.data$group))

  merge_similar_clusters <- function(profile_df, threshold = 0.9) {
    if (nrow(profile_df) == 0) {
      return(tibble::tibble(
        cluster_id = character(),
        merged_cluster_id = character()
      ))
    }

    mean_tbl <- profile_df |>
      dplyr::group_by(.data$cluster_id, .data$group) |>
      dplyr::summarise(
        v = mean(.data$rlog_expr, na.rm = TRUE),
        .groups = "drop"
      )
    mat_mean <- stats::xtabs(v ~ cluster_id + group, data = mean_tbl)
    ids <- rownames(mat_mean)
    if (length(ids) <= 1) {
      return(tibble::tibble(cluster_id = ids, merged_cluster_id = ids))
    }

    cor_mat <- suppressWarnings(stats::cor(
      t(mat_mean),
      use = "pairwise.complete.obs",
      method = "pearson"
    ))
    cor_mat[is.na(cor_mat)] <- 0

    # Build connected components for pairs with high profile correlation.
    neighbors <- lapply(seq_along(ids), function(i) {
      which(cor_mat[i, ] >= threshold)
    })
    visited <- rep(FALSE, length(ids))
    comp_id <- integer(length(ids))
    comp <- 0L
    for (i in seq_along(ids)) {
      if (visited[i]) {
        next
      }
      comp <- comp + 1L
      q <- i
      visited[i] <- TRUE
      comp_id[i] <- comp
      while (length(q) > 0) {
        cur <- q[[1]]
        q <- q[-1]
        for (nxt in neighbors[[cur]]) {
          if (!visited[nxt]) {
            visited[nxt] <- TRUE
            comp_id[nxt] <- comp
            q <- c(q, nxt)
          }
        }
      }
    }

    pick_rep <- function(v) {
      nums <- suppressWarnings(as.integer(sub(
        "^.*_([0-9]+)$",
        "\\1",
        v,
        perl = TRUE
      )))
      if (all(is.na(nums))) {
        return(sort(v)[[1]])
      }
      v[[order(ifelse(is.na(nums), Inf, nums), v)[[1]]]]
    }

    comp_tbl <- tibble::tibble(cluster_id = ids, comp = comp_id) |>
      dplyr::group_by(.data$comp) |>
      dplyr::mutate(merged_cluster_id = pick_rep(.data$cluster_id)) |>
      dplyr::ungroup() |>
      dplyr::select("cluster_id", "merged_cluster_id")

    comp_tbl
  }

  merge_map <- merge_similar_clusters(profile_long, threshold = 0.9)

  cluster_gene <- cluster_gene |>
    dplyr::left_join(merge_map, by = "cluster_id") |>
    dplyr::mutate(
      cluster_id = dplyr::if_else(
        !is.na(.data$merged_cluster_id) & .data$merged_cluster_id != "",
        .data$merged_cluster_id,
        .data$cluster_id
      )
    ) |>
    dplyr::select(-dplyr::any_of("merged_cluster_id")) |>
    dplyr::group_by(.data$cluster_id, .data$ens_gene) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()

  profile_long <- profile_long |>
    dplyr::left_join(merge_map, by = "cluster_id") |>
    dplyr::mutate(
      cluster_id = dplyr::if_else(
        !is.na(.data$merged_cluster_id) & .data$merged_cluster_id != "",
        .data$merged_cluster_id,
        .data$cluster_id
      )
    ) |>
    dplyr::select(-dplyr::any_of("merged_cluster_id"))

  cluster_summary <- cluster_gene |>
    dplyr::group_by(.data$cluster_id) |>
    dplyr::summarise(
      gene_count = dplyr::n_distinct(.data$ens_gene),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$cluster_id)

  group_levels <- unique(as.character(sample_meta$group))
  cluster_profile <- profile_long |>
    dplyr::mutate(group = factor(.data$group, levels = group_levels))

  # 注: 図 (ggplot) はレポートが DEGCluster_profile から描くので、ここでは作らない。
  list(
    cluster_gene = cluster_gene,
    cluster_summary = cluster_summary,
    cluster_profile = cluster_profile,
    plot = NULL,
    merge_map = merge_map
  )
}
