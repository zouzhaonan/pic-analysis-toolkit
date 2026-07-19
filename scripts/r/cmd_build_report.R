#!/usr/bin/env Rscript

# 役割:
#   解析出力ディレクトリを読み取り、自己完結型の HTML レポートを生成する。
# 入力:
#   --out-dir <dir>      : mapping / deseq2 / enrich の出力ベース (必須)
# 出力:
#   <out-dir>/report_<project>.html (project = genome_run ごとに 1 ファイル)

script_path <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path)))
asset_dir <- normalizePath(file.path(script_dir, "..", "assets"), winslash = "/", mustWork = TRUE)

source(file.path(script_dir, "analysis_plot_spec.R"))
source(file.path(script_dir, "plot_common.R"))
source(file.path(script_dir, "enrichment_plots.R"))
source(file.path(script_dir, "report_build.R"))

parse_report_args <- function(args) {
  out <- list(out_dir = NULL)
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--out-dir") { out$out_dir <- args[[i + 1]]; i <- i + 2 }
    else { stop(sprintf("Unknown option: %s", a), call. = FALSE) }
  }
  if (is.null(out$out_dir)) stop("Missing required option: --out-dir", call. = FALSE)
  out
}

main <- function() {
  args <- parse_report_args(commandArgs(trailingOnly = TRUE))
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
  })

  out_dir <- normalizePath(args$out_dir, mustWork = TRUE)

  # xenograft 統合レポート (分類サマリ + graft/host) を優先判定
  if (pic_is_xenograft_out(out_dir)) {
    message("[INFO] report: xenograft integrated (classification + graft/host)")
    html <- build_xenograft_report(out_dir, asset_dir)
    message(sprintf("[INFO] report written: %s", html))
    return(invisible(NULL))
  }

  projects <- pic_report_discover_projects(out_dir)
  if (length(projects) == 0) {
    stop(sprintf("deseq2 の stats_*.csv が見つかりません: %s", out_dir), call. = FALSE)
  }

  # mapping_sum (共有)
  msum <- NULL
  msum_files <- pic_list_summary(out_dir, "^mapping_sum__.*\\.tsv$")
  if (length(msum_files) > 0) {
    msum <- read_mapping_sum(msum_files[[1]])
  }

  for (desc in projects) {
    message(sprintf("[INFO] report: project=%s", desc$project))
    html <- build_report_for_project(desc, out_dir, msum, asset_dir)
    message(sprintf("[INFO] report written: %s", html))
  }
}

main()
