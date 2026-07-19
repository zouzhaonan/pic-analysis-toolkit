# 役割:
#   PCA と sample correlation の元データ (CSV) を出力する。
#   静的な散布図 / 相関ヒートマップ画像はレポート (plotly / HTML) 側で CSV から
#   動的生成するため、ここでは生成しない。
# 入力:
#   DESeq2 object と sample label。
# 出力:
#   PCA_RegLog_<project>.csv, PCA_Variance_<project>.csv, correlation_<project>.csv

build_sample_correlation <- function(dds) {
  normalized_matrix <- DESeq2::counts(dds, normalized = TRUE)
  stats::cor(normalized_matrix, method = "pearson", use = "pairwise.complete.obs")
}

save_correlation_outputs <- function(correlation_matrix, out_dir, project_name) {
  readr::write_csv(
    tibble::rownames_to_column(as.data.frame(correlation_matrix), var = "sample"),
    file.path(out_dir, sprintf("correlation_%s.csv", project_name))
  )
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
  pca_plot_dir <- file.path(out_dir, sprintf("PCA_%s", project_name))
  dir.create(pca_plot_dir, recursive = TRUE, showWarnings = FALSE)
  pca_df <- tibble::as_tibble(pca_res$x, rownames = "sample") |>
    dplyr::left_join(label_tbl, by = "sample")

  readr::write_csv(pca_df, file.path(pca_plot_dir, sprintf("PCA_RegLog_%s.csv", project_name)))
  readr::write_csv(tibble::rownames_to_column(pca_var, "metric"), file.path(pca_plot_dir, sprintf("PCA_Variance_%s.csv", project_name)))
  save_correlation_outputs(build_sample_correlation(dds), out_dir, project_name)
}
