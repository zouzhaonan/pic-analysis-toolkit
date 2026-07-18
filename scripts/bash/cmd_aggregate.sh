#!/usr/bin/env bash

# 役割:
#   遺伝子全体 (TSS→TES) の aggregation / metagene profile を作成する
#   `pic aggregate`。pic mapping が生成した CPM 正規化 bigwig
#   (bw/<count_prefix>_umi.cpm.bw) を再利用し、UCSC bigWigAverageOverBed で
#   各遺伝子を N 分割したビンの平均を高速取得、R で全遺伝子平均のメタジーン
#   プロファイルを計算する (deepTools 非依存)。
# 入力:
#   pic mapping/all の出力ディレクトリ (bw/ と deftable_<run>_<genome>.tsv)。
# 出力:
#   <out>/aggregate_profile_<run>_<genome>.csv  (sample,group,binpos,pos,value; プロジェクト直下)
# 注記:
#   実体は R (scripts/r/cmd_aggregate.R)。bigWigAverageOverBed は pic の conda
#   パッケージに同梱される前提。GTF は ${PIC_LIB} 配下を参照する。

run_aggregate_subcommand() {
  local arg
  for arg in "$@"; do
    if [[ "$arg" == "--help" ]]; then
      show_help_for_subcommand aggregate
      return 0
    fi
  done
  if [[ $# -eq 0 ]]; then
    show_help_for_subcommand aggregate
    return 0
  fi

  debug_log "STEP=run_aggregate_subcommand"
  debug_kv "input.rscript" "${PIC_R_DIR}/cmd_aggregate.R"
  debug_kv "input.args" "$*"
  Rscript "${PIC_R_DIR}/cmd_aggregate.R" "$@"
}
