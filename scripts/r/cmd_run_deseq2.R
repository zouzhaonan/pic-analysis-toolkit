#!/usr/bin/env Rscript

# 役割:
#   genome ごとの deftable から DESeq2 -> plots を実行する。
# 入力:
#   deftable、genome 名、出力先、各種解析パラメータ。
# 出力:
#   DESeq2 結果、各種集計表、plots。

script_path <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path)))
analysis_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "analysis_plot_spec.R"))
source(file.path(script_dir, "analysis_runtime.R"))
paths <- pic_runtime_paths(analysis_root)
biomart_lookup_file <- paths$biomart_lookup_file

source(file.path(script_dir, "deseq2_common.R"))
source(file.path(script_dir, "biomart_cache.R"))
source(file.path(script_dir, "deseq2_data.R"))
source(file.path(script_dir, "plot_common.R"))
source(file.path(script_dir, "contrast_plots.R"))
source(file.path(script_dir, "dimension_plots.R"))

# --part (help/help_deseq2.txt):
#   all     … 全部 (既定)。threads > 1 なら rlog を fork した子で DESeq() と同時に回す。
#   fit     … カウント → サイズ因子 → DESeq → results → 表 (stats・UMI_count・normalizedCountTable・
#              DEG_count・DEG_normalizedCountTable・DEGList) と analysis_params。
#   rlog    … カウント → サイズ因子 → rlog(blind) → PCA・相関。fit を待たずに (同時に) 回せる。
#   cluster … fit と rlog の中間ファイルから DEGCluster_*。
# 部分の済んだ印と中間ファイルは <out-dir>/.pic_deseq2_parts_<project>/ (deseq2_common.R)。

# fit の表を書く。書いたファイルのパスと stats を返す。
write_fit_outputs <- function(args, deg_dir, mat, fit, contrasts, e2g, label) {
  # 採用したサイズ因子の方法をレポート用に記録する。
  # サンプルシートに複数 genome があると本関数は genome ごとに呼ばれるため、
  # キーに genome を含めないと後の genome が前の記録を上書きしてしまう。
  pic_write_analysis_param(
    args$out_dir, paste0("size_factor_method__", args$genome),
    if (exists("pic_size_factor_method", envir = globalenv())) {
      get("pic_size_factor_method", envir = globalenv())
    } else "poscounts"
  )
  pic_write_analysis_param(
    args$out_dir, paste0("dispersion_method__", args$genome),
    if (exists("pic_dispersion_method", envir = globalenv())) {
      get("pic_dispersion_method", envir = globalenv())
    } else "group"
  )
  message("[INFO] Building stats and summary tables")
  stats_tables <- build_stats_tables(fit$stats, contrasts, e2g, args$fdr)
  normcount <- build_normalized_counts(fit$dds, e2g)
  deg_count_summary <- build_deg_count_summary(stats_tables$stats, contrasts, args$fdr)

  # 注: deftable のコピー / Num_UMIs_genes は下流・レポートのいずれからも参照されないため書き出さない
  #     (ルートの deftable と mapping_sum で代替される)。
  files <- c(
    file.path(args$out_dir, sprintf("UMI_count_%s.csv", args$project_name)),
    file.path(args$out_dir, sprintf("stats_%s.csv", args$project_name)),
    file.path(args$out_dir, sprintf("normalizedCountTable_%s.csv", args$project_name))
  )
  pic_write_csv_atomic(mat, files[[1]])
  pic_write_csv_atomic(stats_tables$stats, files[[2]])
  pic_write_csv_atomic(normcount, files[[3]])

  if (nrow(deg_count_summary) > 0) {
    f <- file.path(deg_dir, sprintf("DEG_count_%s.csv", args$project_name))
    pic_write_csv_atomic(deg_count_summary, f)
    files <- c(files, f)
  }

  message("[INFO] Saving contrast-level tables")
  files <- c(files, save_contrast_plots(
    dds = fit$dds,
    label = label,
    stats = stats_tables$stats,
    contrasts = contrasts,
    normalized_count_table = normcount,
    out_dir = deg_dir,
    contrast_plot_dir = args$out_dir,
    project_name = args$project_name,
    fdr = args$fdr
  ))
  list(stats = stats_tables$stats, files = files)
}

# DEGCluster_* を書く。書いたファイルのパスを返す。
write_cluster_outputs <- function(args, deg_dir, degpattern) {
  if (nrow(degpattern$cluster_gene) > 0) {
    degpattern_dir <- file.path(deg_dir, sprintf("DEGCluster_%s", args$project_name))
    dir.create(degpattern_dir, recursive = TRUE, showWarnings = FALSE)
    files <- file.path(degpattern_dir, sprintf(c("DEGCluster_gene_for_ora_%s.csv", "DEGCluster_summary_%s.csv",
                                                 "DEGCluster_profile_%s.csv"), args$project_name))
    pic_write_csv_atomic(degpattern$cluster_gene, files[[1]])
    pic_write_csv_atomic(degpattern$cluster_summary, files[[2]])
    pic_write_csv_atomic(degpattern$cluster_profile, files[[3]])
    # 注: DEGCluster_merge_map はどこからも参照されないため書き出さない
    return(files)
  }
  message("[INFO] DEG clustering returned no cluster assignments")
  character(0)
}

save_rds_atomic <- function(x, path) {
  tmp <- pic_tmp_path(path)
  saveRDS(x, tmp, compress = FALSE)
  pic_commit_tmp(tmp, path)
}

load_count_matrix <- function(args, def, e2g) {
  message("[INFO] Building count matrix")
  umi <- load_umi_tables(def, args$count_dir)
  build_count_matrix(umi, e2g)
}

run_part_fit <- function(args, deg_dir, def, contrasts, e2g, key, parts_dir) {
  begin_part(parts_dir, "fit")
  mat <- load_count_matrix(args, def, e2g)
  prep <- build_deseq_dds(mat, def)
  message("[INFO] Running DESeq2")
  fit <- fit_deseq(prep$dds, contrasts, threads = args$threads, mem_gb = args$mem_gb)
  out <- write_fit_outputs(args, deg_dir, mat, fit, contrasts, e2g, prep$label)
  f_stats <- file.path(parts_dir, "fit_stats.rds")
  save_rds_atomic(as.data.frame(out$stats[, c("ens_gene", "padj")]), f_stats)
  write_part_marker(parts_dir, "fit", key, c(out$files, f_stats), fdr = args$fdr)
  message(sprintf("[INFO] --part fit completed (%s)", part_marker_path(parts_dir, "fit")))
}

run_part_rlog <- function(args, def, e2g, key, parts_dir) {
  begin_part(parts_dir, "rlog")
  mat <- load_count_matrix(args, def, e2g)
  prep <- build_deseq_dds(mat, def)
  message("[INFO] Computing rlog (blind)")
  rlog_mat <- compute_rlog_matrix(prep$dds)
  if (is.null(rlog_mat)) rlog_mat <- SummarizedExperiment::assay(DESeq2::rlog(prep$dds))  # 失敗ならエラーで止まる (all と同じ)
  message("[INFO] Running PCA/correlation")
  files <- run_dimension_reduction(
    dds = prep$dds, label = prep$label, out_dir = args$out_dir,
    project_name = args$project_name, ntop = args$ntop, rlog_mat = rlog_mat, e2g = e2g
  )
  f_rlog <- file.path(parts_dir, "rlog.rds")
  save_rds_atomic(list(rlog_mat = rlog_mat, label = prep$label), f_rlog)
  write_part_marker(parts_dir, "rlog", key, c(files, f_rlog))
  message(sprintf("[INFO] --part rlog completed (%s)", part_marker_path(parts_dir, "rlog")))
}

run_part_cluster <- function(args, deg_dir, e2g, key, parts_dir) {
  require_part(parts_dir, "fit", key, args, fdr = args$fdr)
  require_part(parts_dir, "rlog", key, args)
  begin_part(parts_dir, "cluster")
  stats_min <- readRDS(file.path(parts_dir, "fit_stats.rds"))
  rl <- readRDS(file.path(parts_dir, "rlog.rds"))
  message("[INFO] Running DEG clustering on rlog matrix")
  degpattern <- build_degpattern_outputs(
    dds = NULL, label = rl$label, stats = stats_min, e2g = e2g, fdr = args$fdr, rlog_mat = rl$rlog_mat
  )
  files <- write_cluster_outputs(args, deg_dir, degpattern)
  write_part_marker(parts_dir, "cluster", key, c(files, part_marker_path(parts_dir, c("fit", "rlog"))), fdr = args$fdr)
  message(sprintf("[INFO] --part cluster completed (%s)", part_marker_path(parts_dir, "cluster")))
}

run_part_all <- function(args, deg_dir, def, contrasts, e2g, parts_dir) {
  # 部分ごとの印は、all の結果と食い違わないように消す (all は中間ファイルを作らない)。
  unlink(parts_dir, recursive = TRUE)
  mat <- load_count_matrix(args, def, e2g)
  prep <- build_deseq_dds(mat, def)
  mem <- deseq2_mem_estimate_gb(nrow(prep$dds), ncol(prep$dds))
  master_pid <- Sys.getpid()
  # rlog の子は worker 1 本 (親のほかに 1 スレッドと rlog の分のメモリ) として数える。
  rlog_job <- if (.Platform$OS.type == "unix" && args$threads >= 2L &&
                  args$mem_gb - mem$main >= mem$rlog_worker) {
    start_rlog_job(prep$dds, master_pid)
  } else NULL
  message("[INFO] Running DESeq2")
  fit <- fit_deseq(prep$dds, contrasts, threads = args$threads, mem_gb = args$mem_gb, master_pid = master_pid,
                   busy_threads = if (is.null(rlog_job)) 0L else 1L,
                   busy_gb = if (is.null(rlog_job)) 0 else mem$rlog_worker)
  out <- write_fit_outputs(args, deg_dir, mat, fit, contrasts, e2g, prep$label)
  message("[INFO] Collecting rlog (blind)")
  rlog_mat <- collect_rlog_matrix(fit$dds, rlog_job)
  message("[INFO] Running DEG clustering on rlog matrix")
  degpattern <- build_degpattern_outputs(
    dds = fit$dds, label = prep$label, stats = out$stats, e2g = e2g, fdr = args$fdr, rlog_mat = rlog_mat
  )
  write_cluster_outputs(args, deg_dir, degpattern)
  message("[INFO] Running PCA/correlation")
  run_dimension_reduction(
    dds = fit$dds, label = prep$label, out_dir = args$out_dir,
    project_name = args$project_name, ntop = args$ntop, rlog_mat = rlog_mat, e2g = e2g
  )
}

main <- function() {
  args <- parse_deseq2_args(commandArgs(trailingOnly = TRUE))
  message(sprintf("[INFO] Loading required packages (--part %s)", args$part))
  load_required_packages(args$part)
  pic_set_tmp_tag(args$part)
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
  # 平坦化 (deseq2/<genome>/ 廃止) に伴い、サブフォルダ名にも <run>_<genome>
  # suffix を付けて複数 genome 混在時の衝突を防ぐ。
  deg_dir <- file.path(args$out_dir, sprintf("DEG_%s", args$project_name))
  dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
  parts_dir <- deseq2_parts_dir(args)
  pic_clean_stale_tmp(c(args$out_dir, file.path(dirname(args$out_dir), "summary")), args$part)

  message("[INFO] Reading deftable and preparing contrasts")
  def <- read_deftable(args$deftable)
  contrasts <- build_contrasts_from_deftable(def)
  key <- deseq2_input_key(args, def)
  message(sprintf("[INFO] Loading biomart cache for genome: %s", args$genome))
  e2g <- load_e2g_from_cache(args$genome, biomart_lookup_file)

  switch(args$part,
    fit = run_part_fit(args, deg_dir, def, contrasts, e2g, key, parts_dir),
    rlog = run_part_rlog(args, def, e2g, key, parts_dir),
    cluster = run_part_cluster(args, deg_dir, e2g, key, parts_dir),
    all = run_part_all(args, deg_dir, def, contrasts, e2g, parts_dir)
  )
  message("[INFO] DESeq2 workflow completed")
}

main()
