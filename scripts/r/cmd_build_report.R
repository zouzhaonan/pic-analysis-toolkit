#!/usr/bin/env Rscript

# 役割:
#   解析出力ディレクトリを読み取り、自己完結型の HTML レポートを生成する。
# 入力:
#   --out-dir <dir>      : mapping / deseq2 / enrich の出力ベース (必須)
#   --rin-run-dir <dir>  : この run の TapeStation の泳動フォルダ (任意・繰り返し可)。渡したときだけ
#                          RNA integrity を出す。フォルダを探したり名前から推測したりはしない
#   --rin-pc-dir <dir>   : RNA ScreenTape 用ポジティブコントロールの泳動フォルダ (任意)
#   --rin-pc-hs-dir <dir>: HS RNA ScreenTape 用ポジティブコントロールの泳動フォルダ (任意)。
#                          最初のサンプル (ladder を除く) を PC に使う
#   --threads <n>        : 並べる worker の数 (任意。既定 1)。節の組み立て・GSEA/ORA の読み・埋め込みファイルの
#                          gzip を fork で並べる。出力は --threads によらず同じ
#   --mem-gb <gb>        : 並べる worker のメモリの上限 (任意。既定は環境変数 PIC_MEM_GB、無ければ制限なし)
# 出力:
#   <out-dir>/report_<project>.html (project = <run>_<genome> ごとに 1 ファイル)

script_path <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path)))
asset_dir <- normalizePath(file.path(script_dir, "..", "assets"), winslash = "/", mustWork = TRUE)

# 文字の扱いを UTF-8 にする (source より前)。LANG が空 / C のままだと、R のコードの「−」「↑」などや
# plotly.min.js の U+2212 (負の目盛りの記号) が enc2utf8 で "<e2><88><92>" という文字列に化ける。
if (!isTRUE(l10n_info()[["UTF-8"]])) {
  for (loc in c("C.UTF-8", "en_US.UTF-8", "ja_JP.UTF-8", "C.utf8", "en_US.utf8")) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc))) && isTRUE(l10n_info()[["UTF-8"]])) break
  }
  if (!isTRUE(l10n_info()[["UTF-8"]])) warning("UTF-8 locale is not available; non-ASCII text in the report may be garbled", call. = FALSE)
}

source(file.path(script_dir, "analysis_plot_spec.R"))
source(file.path(script_dir, "plot_common.R"))
source(file.path(script_dir, "enrichment_plots.R"))
source(file.path(script_dir, "report_build.R"))

parse_report_args <- function(args) {
  out <- list(out_dir = NULL, rin_run_dirs = character(0), rin_pc_dir = NULL, rin_pc_hs_dir = NULL, threads = 1L, mem_gb = NULL)
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--out-dir") { out$out_dir <- args[[i + 1]]; i <- i + 2 }
    else if (a == "--rin-run-dir") { out$rin_run_dirs <- c(out$rin_run_dirs, args[[i + 1]]); i <- i + 2 }
    else if (a == "--rin-dir") {
      stop("--rin-dir is no longer supported: pass each TapeStation run folder with --rin-run-dir <dir> (repeatable) ",
           "and the positive control with --rin-pc-dir / --rin-pc-hs-dir.", call. = FALSE)
    }
    else if (a == "--rin-pc-dir") { out$rin_pc_dir <- args[[i + 1]]; i <- i + 2 }
    else if (a == "--rin-pc-hs-dir") { out$rin_pc_hs_dir <- args[[i + 1]]; i <- i + 2 }
    else if (a == "--threads") { out$threads <- suppressWarnings(as.integer(args[[i + 1]])); i <- i + 2 }
    else if (a == "--mem-gb") { out$mem_gb <- suppressWarnings(as.numeric(args[[i + 1]])); i <- i + 2 }
    else { stop(sprintf("Unknown option: %s", a), call. = FALSE) }
  }
  if (is.null(out$out_dir)) stop("Missing required option: --out-dir", call. = FALSE)
  for (opt in c("--rin-run-dir", "--rin-pc-dir", "--rin-pc-hs-dir")) {
    v <- switch(opt, "--rin-run-dir" = out$rin_run_dirs, "--rin-pc-dir" = out$rin_pc_dir, "--rin-pc-hs-dir" = out$rin_pc_hs_dir)
    for (p in v) if (is.na(p) || !nzchar(p) || !dir.exists(path.expand(p)))
      stop(sprintf("%s: folder not found: %s", opt, p), call. = FALSE)
  }
  if (length(out$threads) != 1 || is.na(out$threads) || out$threads < 1L)
    stop("--threads must be an integer >= 1.", call. = FALSE)
  if (!is.null(out$mem_gb) && (length(out$mem_gb) != 1 || is.na(out$mem_gb) || out$mem_gb <= 0))
    stop("--mem-gb must be a number > 0.", call. = FALSE)
  out
}

main <- function() {
  args <- parse_report_args(commandArgs(trailingOnly = TRUE))
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
  })

  out_dir <- normalizePath(args$out_dir, mustWork = TRUE)
  options(pic.report.rin_run_dirs = path.expand(args$rin_run_dirs))
  if (!is.null(args$rin_pc_dir)) options(pic.report.rin_pc_dir = path.expand(args$rin_pc_dir))
  if (!is.null(args$rin_pc_hs_dir)) options(pic.report.rin_pc_hs_dir = path.expand(args$rin_pc_hs_dir))
  if (any(nzchar(Sys.getenv(c("PIC_RIN_DIR", "PIC_RIN_PC_DIR", "PIC_RIN_PC_HS_DIR")))))
    message("[INFO] RIN: PIC_RIN_DIR / PIC_RIN_PC_DIR / PIC_RIN_PC_HS_DIR are no longer read; ",
            "pass --rin-run-dir / --rin-pc-dir / --rin-pc-hs-dir instead.")
  options(pic.report.threads = args$threads, pic.report.mem_gb = args$mem_gb)

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
