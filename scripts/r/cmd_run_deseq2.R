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

main <- function() {
  args <- parse_deseq2_args(commandArgs(trailingOnly = TRUE))
  message("[INFO] Loading required packages")
  load_required_packages()
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
  deg_dir <- file.path(args$out_dir, "DEG")
  dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

  message("[INFO] Reading deftable and preparing contrasts")
  def <- read_deftable(args$deftable)
  contrasts <- build_contrasts_from_deftable(def)
  message(sprintf("[INFO] Loading biomart cache for genome: %s", args$genome))
  e2g <- load_e2g_from_cache(args$genome, biomart_lookup_file)
  message("[INFO] Building count matrix")
  umi <- load_umi_tables(def, args$count_dir)
  mat <- build_count_matrix(umi, e2g)
  message("[INFO] Running DESeq2")
  deseq_result <- run_deseq(mat, def, contrasts)
  message("[INFO] Building stats and summary tables")
  stats_tables <- build_stats_tables(deseq_result$stats, contrasts, e2g, args$fdr)
  normcount <- build_normalized_counts(deseq_result$dds, e2g)
  num_umi_gene <- build_num_umi_gene(mat)
  deg_count_summary <- build_deg_count_summary(stats_tables$stats, contrasts, args$fdr)
  message("[INFO] Running DEG clustering on rlog matrix")
  degpattern <- build_degpattern_outputs(
    dds = deseq_result$dds,
    label = deseq_result$label,
    stats = stats_tables$stats,
    e2g = e2g,
    fdr = args$fdr
  )

  # 注: deftable のコピー / Num_UMIs_genes は下流・レポートのいずれからも参照されないため書き出さない
  #     (ルートの deftable と mapping_sum で代替される)。
  readr::write_csv(mat, file.path(args$out_dir, sprintf("UMI_count_%s.csv", args$project_name)))
  readr::write_csv(stats_tables$stats, file.path(args$out_dir, sprintf("stats_%s.csv", args$project_name)))
  readr::write_csv(normcount, file.path(args$out_dir, sprintf("normalizedCountTable_%s.csv", args$project_name)))

  if (nrow(deg_count_summary) > 0) {
    readr::write_csv(deg_count_summary, file.path(deg_dir, sprintf("DEG_count_%s.csv", args$project_name)))
  }

  if (nrow(degpattern$cluster_gene) > 0) {
    degpattern_dir <- file.path(deg_dir, "DEGCluster")
    dir.create(degpattern_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(
      degpattern$cluster_gene,
      file.path(degpattern_dir, sprintf("DEGCluster_gene_for_ora_%s.csv", args$project_name))
    )
    readr::write_csv(
      degpattern$cluster_summary,
      file.path(degpattern_dir, sprintf("DEGCluster_summary_%s.csv", args$project_name))
    )
    readr::write_csv(
      degpattern$cluster_profile,
      file.path(degpattern_dir, sprintf("DEGCluster_profile_%s.csv", args$project_name))
    )
    # 注: DEGCluster_merge_map はどこからも参照されないため書き出さない
  } else {
    message("[INFO] DEG clustering returned no cluster assignments")
  }

  message("[INFO] Saving contrast-level plots")
  save_contrast_plots(
    dds = deseq_result$dds,
    label = deseq_result$label,
    stats = stats_tables$stats,
    contrasts = contrasts,
    normalized_count_table = normcount,
    out_dir = deg_dir,
    contrast_plot_dir = args$out_dir,
    project_name = args$project_name,
    fdr = args$fdr
  )

  message("[INFO] Running PCA/correlation plots")
  run_dimension_reduction(
    dds = deseq_result$dds,
    label = deseq_result$label,
    out_dir = args$out_dir,
    project_name = args$project_name,
    ntop = args$ntop
  )

  message("[INFO] DESeq2 workflow completed")
}

main()
