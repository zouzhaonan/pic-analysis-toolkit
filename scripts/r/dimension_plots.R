# 役割:
#   PCA と sample correlation の元データ (CSV) を出力する。
#   静的な散布図 / 相関ヒートマップ画像はレポート (plotly / HTML) 側で CSV から
#   動的生成するため、ここでは生成しない。
# 入力:
#   DESeq2 object と sample label。
# 出力:
#   PCA_RegLog_<project>.csv, PCA_Variance_<project>.csv, PCA_rlogTop_<project>.csv (gene, gene_name, サンプル…),
#   correlation_<project>.csv

build_sample_correlation <- function(dds) {
  normalized_matrix <- DESeq2::counts(dds, normalized = TRUE)
  stats::cor(normalized_matrix, method = "pearson", use = "pairwise.complete.obs")
}

save_correlation_outputs <- function(correlation_matrix, out_dir, project_name) {
  f <- file.path(out_dir, sprintf("correlation_%s.csv", project_name))
  pic_write_csv_atomic(
    tibble::rownames_to_column(as.data.frame(correlation_matrix), var = "sample"),
    f
  )
  invisible(f)
}

# レポートが「群の部分集合だけで PCA を引き直す」ために埋め込む rlog 行列の遺伝子数。
# ntop (=500) より多めに持たせる: 部分集合では分散上位に来る遺伝子が入れ替わるため、
# ブラウザ側で改めて上位 ntop を選び直せる余地が要る。
PCA_POOL_GENES <- 2000L

# 書いたファイルのパスを返す。
# dds は DESeq() の前 (--part rlog。サイズ因子だけ設定したもの) でも後 (--part all) でもよい:
# 相関に使う counts(dds, normalized = TRUE) は カウント / サイズ因子 で、DESeq() はサイズ因子が
# 設定済みなら変えない (README の確かめ)。
run_dimension_reduction <- function(dds, label, out_dir, project_name, ntop, rlog_mat = NULL, e2g = NULL) {
  # rlog_mat は main で 1 回だけ計算した rlog(dds, blind = TRUE) の行列 (rlog の既定も blind = TRUE)。
  # 無ければ (rlog が失敗したときなど) ここで計算し直し、失敗ならそのまま止める (前と同じ)。
  if (is.null(rlog_mat)) rlog_mat <- SummarizedExperiment::assay(DESeq2::rlog(dds))
  rld_df <- as.data.frame(rlog_mat)
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

  f_pca <- file.path(pca_plot_dir, sprintf("PCA_RegLog_%s.csv", project_name))
  f_var <- file.path(pca_plot_dir, sprintf("PCA_Variance_%s.csv", project_name))
  pic_write_csv_atomic(pca_df, f_pca)
  pic_write_csv_atomic(tibble::rownames_to_column(pca_var, "metric"), f_var)

  # 部分集合 PCA 用の rlog 行列。レポートはこれを埋め込み、選択された群のサンプルだけで
  # 分散上位 ntop を選び直し -> 中心化/標準化 -> PCA を引き直す (R の prcomp と同じ手順)。
  pool <- min(PCA_POOL_GENES, nrow(rld_df))
  pool_idx <- rank(-gene_var, ties.method = "first") <= pool
  rld_pool <- round(rld_df[pool_idx, , drop = FALSE], 3)
  f_top <- file.path(pca_plot_dir, sprintf("PCA_rlogTop_%s.csv", project_name))
  # gene (ens_gene) のすぐ右に gene_name (deseq2_data.R の pic_gene_name_of)。
  top_tbl <- tibble::rownames_to_column(as.data.frame(rld_pool), "gene")
  top_tbl <- tibble::add_column(top_tbl, gene_name = pic_gene_name_of(top_tbl$gene, e2g), .after = "gene")
  pic_write_csv_atomic(top_tbl, f_top)

  f_cor <- save_correlation_outputs(build_sample_correlation(dds), out_dir, project_name)
  invisible(c(f_pca, f_var, f_top, f_cor))
}
