# 役割:
#   各 plot script で共通に使う補助関数をまとめる。
# 入力:
#   plot 用 data frame や contrast 情報。
# 出力:
#   plot 作成に使う ggplot object や補助ラベル。

build_plot_theme <- function() {
  palette_name <- pic_plot_spec()$plot$palette_d3

  list(
    ggplot2::theme_minimal(),
    ggsci::scale_color_d3(palette_name)
  )
}

save_plot_png <- function(plot_obj, output_path, width = 10, height = 8) {
  ggplot2::ggsave(output_path, plot = plot_obj, width = width, height = height, dpi = 300, limitsize = FALSE)
}

format_contrast_file_label <- function(display_label) {
  gsub(" / ", "_vs_", display_label, fixed = TRUE)
}

build_contrast_direction_labels <- function(numerator, denominator) {
  c(
    numerator_bias = paste0(numerator, " bias"),
    denominator_bias = paste0(denominator, " bias"),
    not_significant = "not significant"
  )
}

annotate_contrast_points <- function(contrast_stats) {
  contrast_stats |>
    dplyr::mutate(
      gene_label = dplyr::if_else(!is.na(.data$ext_gene) & .data$ext_gene != "", .data$ext_gene, .data$ens_gene)
    )
}

annotate_ma_points <- function(contrast_stats, numerator, denominator) {
  direction_labels <- build_contrast_direction_labels(numerator, denominator)
  ma_score <- with(contrast_stats, abs(log2FoldChange) * log10(baseMean + 1))
  ma_score[!is.finite(ma_score)] <- NA_real_
  score_threshold <- stats::quantile(ma_score, probs = 0.75, na.rm = TRUE, names = FALSE)
  if (!is.finite(score_threshold)) {
    score_threshold <- Inf
  }

  annotate_contrast_points(contrast_stats) |>
    dplyr::mutate(
      ma_score = ma_score,
      direction = dplyr::case_when(
        !is.na(.data$ma_score) & .data$ma_score >= score_threshold & .data$log2FoldChange > 0 ~ direction_labels[["numerator_bias"]],
        !is.na(.data$ma_score) & .data$ma_score >= score_threshold & .data$log2FoldChange < 0 ~ direction_labels[["denominator_bias"]],
        TRUE ~ direction_labels[["not_significant"]]
      )
    )
}

annotate_volcano_points <- function(contrast_stats, numerator, denominator, fdr) {
  direction_labels <- build_contrast_direction_labels(numerator, denominator)

  annotate_contrast_points(contrast_stats) |>
    dplyr::mutate(
      direction = dplyr::case_when(
        !is.na(.data$padj) & .data$padj < fdr & .data$log2FoldChange > 0 ~ direction_labels[["numerator_bias"]],
        !is.na(.data$padj) & .data$padj < fdr & .data$log2FoldChange < 0 ~ direction_labels[["denominator_bias"]],
        TRUE ~ direction_labels[["not_significant"]]
      )
    )
}

select_top_contrast_labels <- function(plot_df, numerator, denominator, ranking_column, decreasing = FALSE) {
  top_n <- pic_plot_spec()$plot$contrast$label_top_n

  direction_labels <- build_contrast_direction_labels(numerator, denominator)

  top_positive <- plot_df |>
    dplyr::filter(.data$direction == direction_labels[["numerator_bias"]], !is.na(.data[[ranking_column]]))
  if (decreasing) {
    top_positive <- dplyr::arrange(top_positive, dplyr::desc(.data[[ranking_column]]), dplyr::desc(.data$log2FoldChange))
  } else {
    top_positive <- dplyr::arrange(top_positive, .data[[ranking_column]], dplyr::desc(.data$log2FoldChange))
  }
  top_positive <- dplyr::slice_head(top_positive, n = top_n)

  top_negative <- plot_df |>
    dplyr::filter(.data$direction == direction_labels[["denominator_bias"]], !is.na(.data[[ranking_column]]))
  if (decreasing) {
    top_negative <- dplyr::arrange(top_negative, dplyr::desc(.data[[ranking_column]]), .data$log2FoldChange)
  } else {
    top_negative <- dplyr::arrange(top_negative, .data[[ranking_column]], .data$log2FoldChange)
  }
  top_negative <- dplyr::slice_head(top_negative, n = top_n)

  dplyr::bind_rows(top_positive, top_negative)
}

add_gene_labels <- function(plot_obj, label_df) {
  if (nrow(label_df) == 0) {
    return(plot_obj)
  }

  plot_obj +
    ggrepel::geom_label_repel(
      data = label_df,
      ggplot2::aes(label = .data$gene_label),
      size = 3,
      max.overlaps = Inf,
      color = "black",
      fill = scales::alpha("white", 0.85),
      label.size = 0.15,
      label.padding = grid::unit(0.12, "lines"),
      segment.color = "#636363",
      show.legend = FALSE
    )
}
