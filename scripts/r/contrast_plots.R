# 役割:
#   contrast ごとの DEG 関連データ (DEG 正規化カウント表, DEG 遺伝子リスト) を出力する。
#   MA / volcano / DEG count / DEG heatmap の静的画像はレポート (plotly) 側で
#   CSV から動的生成するため、ここでは生成しない。
# 入力:
#   stats table, contrasts, normalized counts。
# 出力:
#   DEG_normalizedCountTable_<project>.csv, DEGList/<contrast>__<group>.txt

# DEG (padj < fdr) の正規化カウント表を書き出す (ヒートマップの元データ)。
save_deg_normalized_count <- function(stats, normalized_count_table, out_dir, project_name, fdr) {
  if (is.null(normalized_count_table) || !("ens_gene" %in% colnames(normalized_count_table))) {
    return(invisible(NULL))
  }
  deg_genes <- stats |>
    dplyr::filter(!is.na(.data$padj), .data$padj < fdr) |>
    dplyr::pull(.data$ens_gene) |>
    unique()
  if (length(deg_genes) == 0) {
    return(invisible(NULL))
  }
  deg_normalized_count <- normalized_count_table |>
    dplyr::filter(.data$ens_gene %in% deg_genes)
  if (nrow(deg_normalized_count) == 0) {
    return(invisible(NULL))
  }
  readr::write_csv(
    deg_normalized_count,
    file.path(out_dir, sprintf("DEG_normalizedCountTable_%s.csv", project_name))
  )
}

save_deg_lists <- function(stats, contrasts, out_dir, fdr) {
  deg_list_dir <- file.path(out_dir, "DEGList")
  dir.create(deg_list_dir, recursive = TRUE, showWarnings = FALSE)
  has_ext_gene <- "ext_gene" %in% colnames(stats)

  purrr::iwalk(contrasts, function(contrast_values, contrast_name) {
    numerator <- as.character(contrast_values[[2]])
    denominator <- as.character(contrast_values[[3]])
    display_label <- format_contrast_label(numerator, denominator)
    file_label <- format_contrast_file_label(display_label)

    contrast_deg <- stats |>
      dplyr::filter(
        .data$aspect == display_label,
        !is.na(.data$padj),
        .data$padj < fdr,
        !is.na(.data$log2FoldChange)
      ) |>
      dplyr::mutate(
        ext_gene = if (.env$has_ext_gene) as.character(.data$ext_gene) else NA_character_,
        ens_gene = as.character(.data$ens_gene),
        gene_label = dplyr::if_else(!is.na(.data$ext_gene) & .data$ext_gene != "" & .data$ext_gene != "NA", .data$ext_gene, .data$ens_gene)
      ) |>
      dplyr::filter(!is.na(.data$gene_label), .data$gene_label != "")

    numerator_genes <- contrast_deg |>
      dplyr::filter(.data$log2FoldChange > 0) |>
      dplyr::arrange(.data$padj, dplyr::desc(.data$log2FoldChange)) |>
      dplyr::distinct(.data$gene_label) |>
      dplyr::pull(.data$gene_label)

    denominator_genes <- contrast_deg |>
      dplyr::filter(.data$log2FoldChange < 0) |>
      dplyr::arrange(.data$padj, .data$log2FoldChange) |>
      dplyr::distinct(.data$gene_label) |>
      dplyr::pull(.data$gene_label)

    readr::write_lines(
      numerator_genes,
      file.path(deg_list_dir, sprintf("%s__%s.txt", file_label, numerator))
    )
    readr::write_lines(
      denominator_genes,
      file.path(deg_list_dir, sprintf("%s__%s.txt", file_label, denominator))
    )
  })
}

# contrast 関連の DEG データ (CSV / リスト) を出力する。
# dds / label / contrast_plot_dir は旧・静的画像出力用の名残で、現在は使用しない
# (呼び出し側の互換のため引数は維持)。
save_contrast_plots <- function(dds, label, stats, contrasts, normalized_count_table, out_dir, project_name, fdr, contrast_plot_dir = out_dir) {
  save_deg_normalized_count(
    stats = stats,
    normalized_count_table = normalized_count_table,
    out_dir = out_dir,
    project_name = project_name,
    fdr = fdr
  )

  save_deg_lists(
    stats = stats,
    contrasts = contrasts,
    out_dir = out_dir,
    fdr = fdr
  )
}
