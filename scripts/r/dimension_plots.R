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

# レポートが「群の部分集合だけで PCA を引き直す」ために埋め込む rlog 行列の遺伝子数。
# ntop (=500) より多めに持たせる: 部分集合では分散上位に来る遺伝子が入れ替わるため、
# ブラウザ側で改めて上位 ntop を選び直せる余地が要る。
PCA_POOL_GENES <- 2000L

run_dimension_reduction <- function(dds, label, out_dir, project_name, ntop) {
  rld <- DESeq2::rlog(dds)
  rld_df <- as.data.frame(SummarizedExperiment::assay(rld))
  label_tbl <- tibble::as_tibble(label) |>
    dplyr::select("sample", "group")

  gene_var <- apply(rld_df, 1, var)
  ntop <- min(ntop, nrow(rld_df))
  top_idx <- rank(-gene_var, ties.method = "first") <= ntop
  rld_top <- rld_df[top_idx, , drop = FALSE]

  pca_res <- prcomp(t(rld_top), scale. = TRUE, center = TRUE)
  pca_var <- as.data.frame(summary(pca_res)$importance)
  pca_plot_dir <- file.path(out_dir, sprintf("PCA_%s", project_name))
  dir.create(pca_plot_dir, recursive = TRUE, showWarnings = FALSE)
  pca_df <- tibble::as_tibble(pca_res$x, rownames = "sample") |>
    dplyr::left_join(label_tbl, by = "sample")

  readr::write_csv(pca_df, file.path(pca_plot_dir, sprintf("PCA_RegLog_%s.csv", project_name)))
  readr::write_csv(tibble::rownames_to_column(pca_var, "metric"), file.path(pca_plot_dir, sprintf("PCA_Variance_%s.csv", project_name)))

  # 部分集合 PCA 用の rlog 行列。レポートはこれを埋め込み、選択された群のサンプルだけで
  # 分散上位 ntop を選び直し -> 中心化/標準化 -> PCA を引き直す (R の prcomp と同じ手順)。
  pool <- min(PCA_POOL_GENES, nrow(rld_df))
  pool_idx <- rank(-gene_var, ties.method = "first") <= pool
  rld_pool <- round(rld_df[pool_idx, , drop = FALSE], 3)
  readr::write_csv(tibble::rownames_to_column(as.data.frame(rld_pool), "gene"),
                   file.path(pca_plot_dir, sprintf("PCA_rlogTop_%s.csv", project_name)))

  save_correlation_outputs(build_sample_correlation(dds), out_dir, project_name)
}
