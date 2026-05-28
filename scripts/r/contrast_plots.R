# 役割:
#   contrast ごとの MA, volcano, DEG count, DEG heatmap を作る。
# 入力:
#   stats table, contrasts, normalized counts, rlog matrix。
# 出力:
#   contrast 可視化 PNG と DEG 関連 table。

save_deg_count_plot <- function(deg_count_summary, out_dir, project_name) {
  ccfg <- pic_plot_spec()$plot$contrast
  palette_name <- pic_plot_spec()$plot$palette_d3

  if (nrow(deg_count_summary) == 0) {
    return(invisible(NULL))
  }

  plot_df <- deg_count_summary |>
    dplyr::mutate(biased_group = factor(.data$biased_group))

  plot_obj <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(y = .data$contrast, x = .data$deg_count, fill = .data$biased_group)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$deg_count),
      position = ggplot2::position_dodge(width = 0.8),
      hjust = -0.2,
      size = 3
    ) +
    ggsci::scale_fill_d3(palette_name) +
    ggplot2::labs(x = "DEG count", y = "Contrast", fill = "Biased group") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))

  save_plot_png(
    plot_obj,
    file.path(out_dir, sprintf("DEG_count_%s.png", project_name)),
    width = pic_plot_spec()$plot$png_width,
    height = max(ccfg$deg_count_height_min, dplyr::n_distinct(plot_df$contrast) * ccfg$deg_count_height_per_contrast)
  )
}

save_ma_plot <- function(contrast_stats, out_dir, project_name, display_label, file_label, numerator, denominator) {
  ccfg <- pic_plot_spec()$plot$contrast
  plot_df <- contrast_stats |>
    dplyr::filter(!is.na(.data$baseMean), .data$baseMean > 0, !is.na(.data$log2FoldChange)) |>
    annotate_ma_points(numerator, denominator)

  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }

  direction_labels <- build_contrast_direction_labels(numerator, denominator)
  label_df <- select_top_contrast_labels(plot_df, numerator, denominator, "ma_score", decreasing = TRUE)
  direction_colors <- stats::setNames(
    ccfg$direction_colors,
    c(
      direction_labels[["numerator_bias"]],
      direction_labels[["denominator_bias"]],
      direction_labels[["not_significant"]]
    )
  )

  plot_obj <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$baseMean, y = .data$log2FoldChange, color = .data$direction)
  ) +
    ggplot2::geom_point(size = ccfg$point_size, alpha = ccfg$point_alpha) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = direction_colors) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = ccfg$guide_line_color) +
    ggplot2::labs(title = display_label, x = "baseMean", y = "log2 fold change", color = "Direction") +
    ggplot2::theme_minimal()

  plot_obj <- add_gene_labels(plot_obj, label_df)

  save_plot_png(
    plot_obj,
    file.path(out_dir, sprintf("MA_%s_%s.png", file_label, project_name)),
    width = pic_plot_spec()$plot$png_width,
    height = pic_plot_spec()$plot$png_height
  )
}

save_volcano_plot <- function(contrast_stats, out_dir, project_name, display_label, file_label, numerator, denominator, fdr) {
  ccfg <- pic_plot_spec()$plot$contrast
  plot_df <- contrast_stats |>
    dplyr::filter(!is.na(.data$log2FoldChange), !is.na(.data$padj), .data$padj > 0) |>
    annotate_volcano_points(numerator, denominator, fdr) |>
    dplyr::mutate(
      neg_log10_padj = -log10(.data$padj)
    )

  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }

  direction_labels <- build_contrast_direction_labels(numerator, denominator)
  label_df <- select_top_contrast_labels(plot_df, numerator, denominator, "padj")
  direction_colors <- stats::setNames(
    ccfg$direction_colors,
    c(
      direction_labels[["numerator_bias"]],
      direction_labels[["denominator_bias"]],
      direction_labels[["not_significant"]]
    )
  )

  plot_obj <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$log2FoldChange, y = .data$neg_log10_padj, color = .data$direction)
  ) +
    ggplot2::geom_point(size = ccfg$point_size, alpha = ccfg$point_alpha) +
    ggplot2::scale_color_manual(values = direction_colors) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = ccfg$guide_line_color) +
    ggplot2::labs(title = display_label, x = "log2 fold change", y = "-log10(padj)", color = "Direction") +
    ggplot2::theme_minimal()

  plot_obj <- add_gene_labels(plot_obj, label_df)

  save_plot_png(
    plot_obj,
    file.path(out_dir, sprintf("volcano_%s_%s.png", file_label, project_name)),
    width = pic_plot_spec()$plot$png_width,
    height = pic_plot_spec()$plot$png_height
  )
}

save_deg_heatmap <- function(rlog_matrix, label, stats, normalized_count_table, out_dir, project_name, fdr) {
  hcfg <- pic_plot_spec()$plot$contrast$heatmap
  deg_genes <- stats |>
    dplyr::filter(!is.na(.data$padj), .data$padj < fdr) |>
    dplyr::pull(.data$ens_gene) |>
    unique()

  if (length(deg_genes) == 0) {
    return(invisible(NULL))
  }

  sample_order <- label$sample
  heatmap_matrix <- rlog_matrix[rownames(rlog_matrix) %in% deg_genes, sample_order, drop = FALSE]
  if (nrow(heatmap_matrix) == 0 || ncol(heatmap_matrix) == 0) {
    return(invisible(NULL))
  }
  row_labels <- rownames(heatmap_matrix)
  if (
    !is.null(normalized_count_table) &&
    "ens_gene" %in% colnames(normalized_count_table) &&
    "ext_gene" %in% colnames(normalized_count_table)
  ) {
    label_tbl <- normalized_count_table |>
      dplyr::transmute(
        ens_gene = as.character(.data$ens_gene),
        ext_gene = as.character(.data$ext_gene)
      ) |>
      dplyr::filter(!is.na(.data$ens_gene), .data$ens_gene != "") |>
      dplyr::distinct(.data$ens_gene, .keep_all = TRUE)
    idx <- match(rownames(heatmap_matrix), label_tbl$ens_gene)
    ext <- label_tbl$ext_gene[idx]
    row_labels <- ifelse(!is.na(ext) & ext != "" & ext != "NA", ext, rownames(heatmap_matrix))
  }

  cluster_rows_flag <- nrow(heatmap_matrix) >= 2

  heatmap_obj <- pheatmap::pheatmap(
    heatmap_matrix,
    scale = "row",
    cluster_rows = cluster_rows_flag,
    cluster_cols = FALSE,
    show_rownames = nrow(heatmap_matrix) <= hcfg$rowname_max,
    fontsize_row = hcfg$fontsize_row,
    fontsize_col = hcfg$fontsize_col,
    clustering_distance_rows = "correlation",
    color = grDevices::colorRampPalette(c(hcfg$color_low, hcfg$color_mid, hcfg$color_high))(hcfg$color_steps),
    labels_row = row_labels,
    legend = FALSE,
    silent = TRUE
  )

  grDevices::png(
    filename = file.path(out_dir, sprintf("DEG_heatmap_%s.png", project_name)),
    width = hcfg$width,
    height = hcfg$height,
    units = "in",
    res = 300
  )
  grid::grid.newpage()
  grid::grid.draw(heatmap_obj$gtable)
  grDevices::dev.off()

  ordered_genes <- rownames(heatmap_matrix)
  if (inherits(heatmap_obj$tree_row, "hclust") && !is.null(heatmap_obj$tree_row$order)) {
    ordered_genes <- rownames(heatmap_matrix)[heatmap_obj$tree_row$order]
  }

  if (!is.null(normalized_count_table) && "ens_gene" %in% colnames(normalized_count_table)) {
    deg_normalized_count <- normalized_count_table |>
      dplyr::filter(.data$ens_gene %in% ordered_genes) |>
      dplyr::mutate(ens_gene = factor(.data$ens_gene, levels = ordered_genes)) |>
      dplyr::arrange(.data$ens_gene) |>
      dplyr::mutate(ens_gene = as.character(.data$ens_gene))

    readr::write_csv(
      deg_normalized_count,
      file.path(out_dir, sprintf("DEG_normalizedCountTable_%s.csv", project_name))
    )
  }
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

save_contrast_plots <- function(dds, label, stats, contrasts, normalized_count_table, out_dir, project_name, fdr, contrast_plot_dir = out_dir) {
  ma_dir <- file.path(contrast_plot_dir, "MA")
  volcano_dir <- file.path(contrast_plot_dir, "volcano")
  dir.create(ma_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

  purrr::iwalk(contrasts, function(contrast_values, contrast_name) {
    display_label <- format_contrast_label(contrast_values[[2]], contrast_values[[3]])
    file_label <- format_contrast_file_label(display_label)
    contrast_stats <- stats |>
      dplyr::filter(.data$aspect == display_label)

    save_ma_plot(contrast_stats, ma_dir, project_name, display_label, file_label, contrast_values[[2]], contrast_values[[3]])
    save_volcano_plot(contrast_stats, volcano_dir, project_name, display_label, file_label, contrast_values[[2]], contrast_values[[3]], fdr)
  })

  save_deg_heatmap(
    rlog_matrix = SummarizedExperiment::assay(DESeq2::rlog(dds)),
    label = label,
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
