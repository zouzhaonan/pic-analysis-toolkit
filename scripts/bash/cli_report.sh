#!/usr/bin/env bash

# 役割:
#   解析出力ディレクトリから HTML レポートを生成する `pic report`。
# 入力:
#   --out-dir <dir>  : mapping / deseq2 / enrich の出力ベース (必須)
#   --threads <n>    : 並べる worker の数 (任意。既定 1)
#   --mem-gb <gb>    : worker のメモリの上限 (任意。既定は PIC_MEM_GB)
# 出力:
#   <out-dir>/report_<project>.html

run_report_subcommand() {
  local arg
  for arg in "$@"; do
    if [[ "$arg" == "--help" ]]; then
      show_help_for_subcommand report
      return 0
    fi
  done
  if [[ $# -eq 0 ]]; then
    show_help_for_subcommand report
    return 0
  fi

  debug_log "STEP=run_report_subcommand"
  debug_kv "input.rscript" "${PIC_R_DIR}/cmd_build_report.R"
  debug_kv "input.args" "$*"
  Rscript "${PIC_R_DIR}/cmd_build_report.R" "$@"
}
