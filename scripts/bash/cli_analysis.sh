#!/usr/bin/env bash

run_analysis_subcommand() {
  local subcommand="$1"
  shift || true

  local arg
  for arg in "$@"; do
    if [[ "$arg" == "--help" ]]; then
      show_help_for_subcommand "$subcommand"
      return 0
    fi
  done
  if [[ $# -eq 0 ]]; then
    show_help_for_subcommand "$subcommand"
    return 0
  fi

  case "$subcommand" in
    deseq2)
      debug_log "STEP=run_analysis_subcommand name=deseq2"
      debug_kv "input.rscript" "${PIC_R_DIR}/cmd_run_deseq2.R"
      debug_kv "input.args" "$*"
      Rscript "${PIC_R_DIR}/cmd_run_deseq2.R" "$@"
      ;;
    enrich)
      debug_log "STEP=run_analysis_subcommand name=enrich"
      debug_kv "input.rscript" "${PIC_R_DIR}/cmd_run_enrich.R"
      debug_kv "input.args" "$*"
      Rscript "${PIC_R_DIR}/cmd_run_enrich.R" "$@"
      ;;
    manage-biomart)
      debug_log "STEP=run_analysis_subcommand name=manage-biomart"
      debug_kv "input.rscript" "${PIC_R_DIR}/cmd_manage_biomart.R"
      debug_kv "input.args" "$*"
      Rscript "${PIC_R_DIR}/cmd_manage_biomart.R" "$@"
      ;;
    plot-expression)
      debug_log "STEP=run_analysis_subcommand name=plot-expression"
      debug_kv "input.rscript" "${PIC_R_DIR}/cmd_plot_expression.R"
      debug_kv "input.args" "$*"
      Rscript "${PIC_R_DIR}/cmd_plot_expression.R" "$@"
      ;;
    chipatlas)
      debug_log "STEP=run_analysis_subcommand name=chipatlas"
      debug_kv "input.script" "${PIC_BASH_DIR}/cmd_chipatlas.sh"
      debug_kv "input.args" "$*"
      run_chipatlas_subcommand "$@"
      ;;
  esac
}
