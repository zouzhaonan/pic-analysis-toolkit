# 役割:
#   PCA と sample correlation の可視化をまとめる。
# 入力:
#   DESeq2 object と sample label。
# 出力:
#   PCA と correlation の CSV / PNG。

format_pc_axis_label <- function(pc_name, variance_table) {
  variance_percent <- NA_real_
  if (pc_name %in% colnames(variance_table) && "Proportion of Variance" %in% rownames(variance_table)) {
    variance_percent <- 100 * variance_table["Proportion of Variance", pc_name]
  }

  if (is.na(variance_percent)) {
    return(pc_name)
  }

  sprintf("%s (%d%%)", pc_name, round(variance_percent))
}

make_point_plot <- function(df, x_col, y_col, x_label = x_col, y_label = y_col) {
  dcfg <- pic_plot_spec()$plot$dimension

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], colour = .data$group, label = .data$sample)
  ) +
    ggplot2::geom_point(size = dcfg$pca_point_size, stroke = 1) +
    ggrepel::geom_text_repel(size = dcfg$pca_label_size, max.overlaps = Inf) +
    ggplot2::labs(x = x_label, y = y_label) +
    build_plot_theme()
}

build_sample_correlation <- function(dds) {
  normalized_matrix <- DESeq2::counts(dds, normalized = TRUE)
  stats::cor(normalized_matrix, method = "pearson", use = "pairwise.complete.obs")
}

save_correlation_outputs <- function(correlation_matrix, out_dir, project_name) {
  dcfg <- pic_plot_spec()$plot$dimension
  sample_order <- colnames(correlation_matrix)

  readr::write_csv(
    tibble::rownames_to_column(as.data.frame(correlation_matrix), var = "sample"),
    file.path(out_dir, sprintf("correlation_%s.csv", project_name))
  )

  corr_df <- as.data.frame(as.table(correlation_matrix), stringsAsFactors = FALSE)
  colnames(corr_df) <- c("sample_x", "sample_y", "correlation")
  corr_df$sample_x <- factor(corr_df$sample_x, levels = sample_order)
  corr_df$sample_y <- factor(corr_df$sample_y, levels = rev(sample_order))

  heatmap_plot <- ggplot2::ggplot(
    corr_df,
    ggplot2::aes(x = .data$sample_x, y = .data$sample_y, fill = .data$correlation)
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$correlation)), size = 3) +
    ggplot2::scale_fill_gradient2(
      low = dcfg$corr_color_low,
      mid = dcfg$corr_color_mid,
      high = dcfg$corr_color_high,
      midpoint = 0
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Pearson r") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  save_plot_png(
    heatmap_plot,
    file.path(out_dir, sprintf("correlation_%s.png", project_name)),
    width = max(dcfg$corr_width_min, ncol(correlation_matrix) * dcfg$corr_width_per_sample),
    height = max(dcfg$corr_height_min, ncol(correlation_matrix) * dcfg$corr_height_per_sample)
  )
}

save_pca_pair_plot <- function(pca_df, variance_table, out_dir, project_name, x_pc, y_pc) {
  if (!all(c(x_pc, y_pc) %in% colnames(pca_df))) {
    return(invisible(NULL))
  }

  file_name <- sprintf("PCA_%s_%s_%s.png", x_pc, y_pc, project_name)
  plot_obj <- make_point_plot(
    pca_df,
    x_pc,
    y_pc,
    x_label = format_pc_axis_label(x_pc, variance_table),
    y_label = format_pc_axis_label(y_pc, variance_table)
  )

  save_plot_png(plot_obj, file.path(out_dir, file_name))
}

run_dimension_reduction <- function(dds, label, out_dir, project_name, ntop) {
  rld <- DESeq2::rlog(dds)
  rld_df <- as.data.frame(SummarizedExperiment::assay(rld))
  label_tbl <- tibble::as_tibble(label) |>
    dplyr::select("sample", "group")

  ntop <- min(ntop, nrow(rld_df))
  top_idx <- rank(-apply(rld_df, 1, var), ties.method = "first") <= ntop
  rld_top <- rld_df[top_idx, , drop = FALSE]

  pca_res <- prcomp(t(rld_top), scale. = TRUE, center = TRUE)
  pca_var <- as.data.frame(summary(pca_res)$importance)
  pca_plot_dir <- file.path(out_dir, "PCA")
  dir.create(pca_plot_dir, recursive = TRUE, showWarnings = FALSE)
  pca_df <- tibble::as_tibble(pca_res$x, rownames = "sample") |>
    dplyr::left_join(label_tbl, by = "sample")

  readr::write_csv(pca_df, file.path(pca_plot_dir, sprintf("PCA_RegLog_%s.csv", project_name)))
  readr::write_csv(tibble::rownames_to_column(pca_var, "metric"), file.path(pca_plot_dir, sprintf("PCA_Variance_%s.csv", project_name)))
  save_correlation_outputs(build_sample_correlation(dds), out_dir, project_name)

  save_pca_pair_plot(pca_df, pca_var, pca_plot_dir, project_name, "PC1", "PC2")
  save_pca_pair_plot(pca_df, pca_var, pca_plot_dir, project_name, "PC2", "PC3")
  save_pca_pair_plot(pca_df, pca_var, pca_plot_dir, project_name, "PC3", "PC4")
}
