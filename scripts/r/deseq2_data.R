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
  sample_mat <- umi |>
    dplyr::transmute(
      ens_gene = as.character(.data$ens_gene),
      sample = as.character(.data$sample),
      count = as.numeric(.data$count)
    ) |>
    dplyr::filter(
      !is.na(.data$ens_gene),
      .data$ens_gene != "",
      !is.na(.data$sample),
      .data$sample != ""
    ) |>
    dplyr::group_by(.data$ens_gene, .data$sample) |>
    dplyr::summarise(
      count = sum(.data$count, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = "sample",
      values_from = "count",
      values_fill = 0
    )

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

  e2g_uniq <- e2g |>
    dplyr::transmute(
      ens_gene = as.character(.data$ens_gene),
      ext_gene = as.character(.data$ext_gene),
      biotype = as.character(.data$biotype),
      chr = as.character(.data$chr)
    ) |>
    dplyr::filter(!is.na(.data$ens_gene), .data$ens_gene != "") |>
    dplyr::group_by(.data$ens_gene) |>
    dplyr::summarise(
      ext_gene = {
        v <- .data$ext_gene[!is.na(.data$ext_gene) & .data$ext_gene != ""]
        if (length(v) == 0) NA_character_ else v[[1]]
      },
      biotype = {
        v <- .data$biotype[!is.na(.data$biotype) & .data$biotype != ""]
        if (length(v) == 0) NA_character_ else v[[1]]
      },
      chr = {
        v <- .data$chr[!is.na(.data$chr) & .data$chr != ""]
        if (length(v) == 0) NA_character_ else v[[1]]
      },
      .groups = "drop"
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
  if (!degenerate && length(sf) > 2) {
    # サイズ因子が深度差をどれだけ追随できているかで判定する。
    # 正常なら log(sizeFactor) は log(ライブラリサイズ) と同程度にばらつく。
    # 退化していると、深度が何倍違ってもサイズ因子はほぼ動かない。
    # (異常サンプル 1 本だけが外れる場合に単純な最大/最小比では検知できないため、
    #  ばらつきの比で判定する。)
    sd_sf <- stats::sd(log(sf))
    sd_lib <- stats::sd(log(lib))
    degenerate <- is.finite(sd_lib) && sd_lib > 0.1 && sd_sf < 0.2 * sd_lib
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
  utils::write.table(rows, f, sep = "\t", quote = FALSE, row.names = FALSE)
}

run_deseq <- function(mat, def, contrasts) {
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

    dds <- DESeq2::estimateSizeFactors(dds)
    normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
    base_mean <- rowMeans(normalized_counts, na.rm = TRUE)

    results_list <- purrr::imap(
      contrasts,
      function(contrast_values, contrast_name) {
        numerator <- contrast_values[[2]]
        denominator <- contrast_values[[3]]

        numerator_samples <- label_df$sample[label_df$group == numerator]
        denominator_samples <- label_df$sample[label_df$group == denominator]

        numerator_mean <- rowMeans(
          normalized_counts[, numerator_samples, drop = FALSE],
          na.rm = TRUE
        )
        denominator_mean <- rowMeans(
          normalized_counts[, denominator_samples, drop = FALSE],
          na.rm = TRUE
        )
        log2_fc <- log2((numerator_mean + 1) / (denominator_mean + 1))

        as.data.frame(
          list(
            baseMean = base_mean,
            log2FoldChange = log2_fc,
            lfcSE = NA_real_,
            stat = NA_real_,
            pvalue = NA_real_,
            padj = NA_real_
          ),
          row.names = rownames(normalized_counts)
        )
      }
    )
  } else {
    dds <- deseq_fit

    results_list <- purrr::imap(
      contrasts,
      ~ DESeq2::results(
        dds,
        contrast = .x,
        lfcThreshold = log2(1),
        independentFiltering = FALSE,
        cooksCutoff = FALSE
      )
    )
  }

  list(
    dds = dds,
    label = label_df,
    mat = mat,
    stats = results_list
  )
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

build_rlog_matrix_for_degpatterns <- function(dds) {
  rlog_try <- tryCatch(
    SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE)),
    error = function(err) NULL
  )
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
  min_cluster_size = 5L
) {
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    stop("Missing required R package: ggbeeswarm", call. = FALSE)
  }

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

  rlog_mat <- build_rlog_matrix_for_degpatterns(dds)
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

  vcfg <- pic_plot_spec()$plot$violin
  palette_name <- pic_plot_spec()$plot$palette_d3

  cluster_plot <- ggplot2::ggplot(
    cluster_profile,
    ggplot2::aes(
      x = .data$group,
      y = .data$rlog_expr,
      color = .data$group,
      fill = .data$group
    )
  ) +
    ggplot2::geom_boxplot(
      width = vcfg$box_width,
      outlier.shape = NA,
      alpha = vcfg$box_alpha,
      linewidth = vcfg$box_linewidth
    ) +
    ggbeeswarm::geom_quasirandom(
      method = "pseudorandom",
      size = 0.2,
      alpha = vcfg$point_alpha
    ) +
    ggplot2::facet_wrap(~cluster_id, scales = "free_y") +
    ggsci::scale_color_d3(palette_name, name = "Group") +
    ggsci::scale_fill_d3(palette_name, name = "Group") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = vcfg$legend_title_size),
      legend.text = ggplot2::element_text(size = vcfg$legend_text_size),
      strip.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = vcfg$box_linewidth
      )
    ) +
    ggplot2::labs(y = "rlog expression")

  list(
    cluster_gene = cluster_gene,
    cluster_summary = cluster_summary,
    cluster_profile = cluster_profile,
    plot = cluster_plot,
    merge_map = merge_map
  )
}
