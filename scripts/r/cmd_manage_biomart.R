#!/usr/bin/env Rscript

# 役割:
#   biomart 準備処理の入口だけを持つ。
# 入力:
#   CLI 引数。
# 出力:
#   実処理は sourced script 側で作る。

script_path <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path)))
analysis_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "analysis_runtime.R"))
paths <- pic_runtime_paths(analysis_root)
biomart_lookup_file <- paths$biomart_lookup_file

source(file.path(script_dir, "biomart_shared.R"))
source(file.path(script_dir, "biomart_prepare.R"))

tryCatch(
  main(),
  error = function(err) {
    message(sprintf("Error: %s", conditionMessage(err)))
    cat("\n")
    print_usage()
    quit(save = "no", status = 1)
  }
)
