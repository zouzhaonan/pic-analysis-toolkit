build_contrasts_from_stats <- function(stats) {
  if (!"aspect" %in% colnames(stats)) {
    stop("stats table must include 'aspect' column.", call. = FALSE)
  }

  aspects <- stats |>
    dplyr::filter(!is.na(.data$aspect), .data$aspect != "") |>
    dplyr::distinct(.data$aspect) |>
    dplyr::pull(.data$aspect)

  parsed <- lapply(aspects, function(aspect_label) {
    parts <- strsplit(aspect_label, " / ", fixed = TRUE)[[1]]
    if (length(parts) != 2) {
      return(NULL)
    }
    c("group", parts[[1]], parts[[2]])
  })

  parsed <- Filter(Negate(is.null), parsed)
  names(parsed) <- vapply(parsed, function(x) paste0(x[[2]], "_vs_", x[[3]]), character(1))

  parsed
}

# GSEA の順位の表 (gene_id = gene_col の値、score、その gene_id に選ばれた行の ens_gene)。
# 同じ gene_id が複数の行 (別の ens_gene) にあるときは |score| の最も大きい 1 行を使う。
ranked_gene_table <- function(contrast_stats, gene_col = "ens_gene") {
  if (!gene_col %in% colnames(contrast_stats)) {
    stop(sprintf("stats table must include '%s' column for GSEA.", gene_col), call. = FALSE)
  }
  ens <- if ("ens_gene" %in% colnames(contrast_stats)) as.character(contrast_stats$ens_gene) else rep(NA_character_, nrow(contrast_stats))

  if (all(c("stat", "log2FoldChange") %in% colnames(contrast_stats))) {
    rank_df <- contrast_stats |>
      dplyr::transmute(
        gene_id = .data[[gene_col]],
        score = dplyr::if_else(
          !is.na(as.numeric(.data$stat)),
          as.numeric(.data$stat),
          as.numeric(.data$log2FoldChange)
        )
      )
  } else if ("stat" %in% colnames(contrast_stats)) {
    rank_df <- contrast_stats |>
      dplyr::transmute(
        gene_id = .data[[gene_col]],
        score = as.numeric(.data$stat)
      )
  } else if ("log2FoldChange" %in% colnames(contrast_stats)) {
    rank_df <- contrast_stats |>
      dplyr::transmute(
        gene_id = .data[[gene_col]],
        score = as.numeric(.data$log2FoldChange)
      )
  } else {
    stop(
      "stats table must include either 'stat' or 'log2FoldChange' for GSEA.",
      call. = FALSE
    )
  }
  rank_df$ens_gene <- ens

  rank_df |>
    dplyr::filter(!is.na(.data$gene_id), .data$gene_id != "", is.finite(.data$score)) |>
    dplyr::group_by(.data$gene_id) |>
    dplyr::slice_max(order_by = abs(.data$score), n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(.data$score))
}

# ranked_gene_table から GSEA に渡す順位のベクトルを作る。
ranked_list_from_table <- function(rank_df) {
  ranked <- stats::setNames(rank_df$score, rank_df$gene_id)
  ranked <- sort(ranked, decreasing = TRUE)

  # clusterProfiler::gseGO may fail on some inputs when many tied scores exist.
  # Add a tiny monotonic offset to preserve order while making values strictly decreasing.
  adjusted_scores <- as.numeric(ranked) - seq_along(ranked) * 1e-12
  stats::setNames(adjusted_scores, names(ranked))
}

build_ranked_gene_list <- function(contrast_stats, gene_col = "ens_gene") {
  ranked_list_from_table(ranked_gene_table(contrast_stats, gene_col))
}

# 1 つのセルに "/" で並んだ遺伝子 (geneID・core_enrichment) を、同じ並びの ID に置き換える。
# id_map: 名前が遺伝子 (symbol)、値が ID (1 つの symbol に複数の ens_gene があれば ";" でつないだもの)。
# 引けない遺伝子は空にする (並びの位置は保つ)。NA のセルは NA のまま。
pic_genes_to_ids <- function(x, id_map) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x)
  if (!any(ok)) return(out)
  out[ok] <- vapply(strsplit(x[ok], "/", fixed = TRUE), function(g) {
    v <- if (is.null(id_map)) rep(NA_character_, length(g)) else unname(id_map[g])
    v[is.na(v)] <- ""
    paste(v, collapse = "/")
  }, character(1))
  out
}

calc_str <- function(n) {
  result <- {
    parts <- do.call(rbind, strsplit(n, "/", fixed = TRUE))
    as.numeric(parts[,1]) / as.numeric(parts[,2])
  }

  result
}

# symbol -> その symbol に当たる ens_gene (複数なら ";" でつなぐ。出てきた順)。名前は unique(symbol) の順。
symbol_to_ids <- function(symbol, ens) {
  u <- unique(symbol)
  sp <- split(as.character(ens), factor(symbol, levels = u))
  stats::setNames(vapply(sp, function(z) paste(unique(z[!is.na(z) & z != ""]), collapse = ";"), character(1)), u)
}

load_precomputed_deg_clusters <- function(cluster_file, stats, gene_col = "human_gene_symbol") {
  if (is.null(cluster_file) || !nzchar(as.character(cluster_file))) {
    return(NULL)
  }
  if (!file.exists(cluster_file)) {
    stop(sprintf("degPattern cluster file was not found: %s", cluster_file), call. = FALSE)
  }

  df <- suppressMessages(
    readr::read_csv(cluster_file, show_col_types = FALSE, progress = FALSE)
  )
  if (nrow(df) == 0) return(tibble::tibble(cluster_id = character(), genes = list(), gene_count = integer(), gene_ids = list()))
  if (!"cluster_id" %in% colnames(df)) {
    stop("degPattern cluster file must include 'cluster_id' column.", call. = FALSE)
  }

  preferred_cols <- unique(c(
    gene_col,
    if (identical(gene_col, "source_gene_symbol")) c("source_gene_symbol", "ext_gene") else character(),
    if (identical(gene_col, "human_gene_symbol")) c("human_gene_symbol", "human_ortholog") else character()
  ))
  preferred_cols <- preferred_cols[preferred_cols %in% colnames(df)]

  # ens_gene は ORA の geneID を ID に戻す (geneID_ens_gene) のに使う。
  df_ens <- if ("ens_gene" %in% colnames(df)) as.character(df$ens_gene) else rep(NA_character_, nrow(df))
  if (length(preferred_cols) > 0) {
    use_col <- preferred_cols[[1]]
    mapped <- df |>
      dplyr::transmute(
        cluster_id = as.character(.data$cluster_id),
        human_gene_symbol = as.character(.data[[use_col]]),
        ens_gene = df_ens
      )
  } else if ("human_gene_symbol" %in% colnames(df)) {
    mapped <- df |>
      dplyr::transmute(
        cluster_id = as.character(.data$cluster_id),
        human_gene_symbol = as.character(.data$human_gene_symbol),
        ens_gene = df_ens
      )
  } else if ("ens_gene" %in% colnames(df)) {
    if (!all(c("ens_gene", gene_col) %in% colnames(stats))) {
      stop(sprintf("stats must include ens_gene and %s to map precomputed clusters.", gene_col), call. = FALSE)
    }
    map_tbl <- stats |>
      dplyr::transmute(ens_gene = as.character(.data$ens_gene), human_gene_symbol = as.character(.data[[gene_col]])) |>
      dplyr::filter(!is.na(.data$ens_gene), .data$ens_gene != "", !is.na(.data$human_gene_symbol), .data$human_gene_symbol != "") |>
      dplyr::distinct(.data$ens_gene, .data$human_gene_symbol)
    mapped <- df |>
      dplyr::transmute(cluster_id = as.character(.data$cluster_id), ens_gene = as.character(.data$ens_gene)) |>
      dplyr::left_join(map_tbl, by = "ens_gene") |>
      dplyr::select("cluster_id", "human_gene_symbol", "ens_gene")
  } else {
    stop(
      sprintf(
        "degPattern cluster file must include one of {%s, ens_gene}.",
        paste(unique(c(gene_col, "ext_gene", "human_gene_symbol")), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  mapped |>
    dplyr::filter(!is.na(.data$cluster_id), .data$cluster_id != "", !is.na(.data$human_gene_symbol), .data$human_gene_symbol != "") |>
    dplyr::group_by(.data$cluster_id) |>
    dplyr::summarise(
      genes = list(unique(.data$human_gene_symbol)),
      gene_count = dplyr::n_distinct(.data$human_gene_symbol),
      gene_ids = list(symbol_to_ids(.data$human_gene_symbol, .data$ens_gene)),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$cluster_id)
}
