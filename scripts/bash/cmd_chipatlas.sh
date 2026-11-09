#!/usr/bin/env bash

# pic chipatlas: ChIP-Atlas の EA / PAGE を DDBJ の WABI に投げる。本体は scripts/r/cmd_chipatlas.R。
# --detach を付けると、裏 (setsid nohup) で走らせて、すぐに戻る。ログは <out-dir>/run/chipatlas.log。

run_chipatlas_subcommand() {
  local -a args=()
  local detach=0
  local out_dir=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
    --help)
      show_help_for_subcommand "chipatlas"
      return 0
      ;;
    --detach)
      detach=1
      shift
      ;;
    --out-dir)
      out_dir="${2:-}"
      args+=("$1" "${2:-}")
      shift 2 || shift
      ;;
    --deg-list-dir)
      handle_error "--deg-list-dir is no longer used. Pass --stats <stats_<project>.csv> (DEGs are taken from it with --fdr)."
      return 1
      ;;
    *)
      args+=("$1")
      shift
      ;;
    esac
  done

  if [[ "$detach" -eq 1 ]]; then
    if [[ -z "$out_dir" ]]; then
      handle_error "--out-dir is required."
      return 1
    fi
    mkdir -p "$out_dir/run"
    local log_file="$out_dir/run/chipatlas.log"
    setsid nohup Rscript "${PIC_R_DIR}/cmd_chipatlas.R" "${args[@]}" >>"$log_file" 2>&1 </dev/null &
    log_info "[ChIP-Atlas] Running in background (pid $!). Log: $log_file"
    return 0
  fi

  Rscript "${PIC_R_DIR}/cmd_chipatlas.R" "${args[@]}"
}
