#!/usr/bin/env bash

# 役割:
#   コマンドごとに、pic の作業フォルダを作り直したり消したりする。
# 入力:
#   init_config で設定された各種グローバル変数。
# 出力:
#   処理ごとに必要な作業フォルダや出力ファイルを作成・削除する。

prepare_workspace_dirs() {
  # log/ は診断ログを実行中から直接書き込む先。tmp/ は scratch 専用。
  mkdir -p "$OUTPUT_DIR" "$TMP_DIR" "$LOG_DIR" "$DEFTABLE_DIR"
}

is_protected_workspace_path() {
  local candidate="$1"
  local c abs_raw abs_demux abs_sheet

  c="$(resolve_existing_path "$candidate" 2>/dev/null || true)"
  [[ -z "$c" ]] && return 1

  abs_raw="$(resolve_existing_path "$RAW_FASTQ_DIR" 2>/dev/null || true)"
  abs_demux="$(resolve_existing_path "$DEMUX_FASTQ_DIR" 2>/dev/null || true)"
  abs_sheet="$(resolve_existing_path "$SAMPLE_SHEET" 2>/dev/null || true)"

  [[ -n "$abs_raw" && "$c" = "$abs_raw" ]] && return 0
  [[ -n "$abs_demux" && "$c" = "$abs_demux" ]] && return 0
  [[ -n "$abs_sheet" && "$c" = "$abs_sheet" ]] && return 0
  return 1
}

safe_rm_rf() {
  local target="$1"
  [[ -z "$target" ]] && return 0
  [[ ! -e "$target" ]] && return 0

  if is_protected_workspace_path "$target"; then
    debug_log "SKIP_REMOVE protected_path=$(resolve_existing_path "$target")"
    return 0
  fi
  rm -fr "$target"
}

reset_demux_outputs() {
  safe_rm_rf "$TMP_DIR"

  if [[ "$SKIP_DEMUX" = 0 ]]; then
    safe_rm_rf "$DEMUX_FASTQ_DIR"
  fi
}

reset_mapping_outputs() {
  safe_rm_rf "$TMP_DIR"
  safe_rm_rf "$BAM_DIR"
  safe_rm_rf "$COUNTS_DIR"
  safe_rm_rf "$LOG_DIR"
  safe_rm_rf "$MAP_SUMMARY"
  safe_rm_rf "${OUTPUT_DIR}/mapping_sum.tsv"
}

reset_bigwig_outputs() {
  safe_rm_rf "$BIGWIG_DIR"
  safe_rm_rf "$UMI_BAM_DIR"
  safe_rm_rf "${OUTPUT_DIR}/sim_umis.png"
}

prepare_workspace_for_run() {
  reset_demux_outputs
  reset_mapping_outputs
  reset_bigwig_outputs
  prepare_workspace_dirs

  mkdir -p "$COUNTS_DIR" "$BAM_DIR" "$BIGWIG_DIR" "$UMI_BAM_DIR"

  if [[ "$SKIP_DEMUX" = 0 ]]; then
    mkdir -p "$DEMUX_FASTQ_DIR"
  fi
}

prepare_workspace_for_demux() {
  reset_demux_outputs
  prepare_workspace_dirs

  if [[ "$SKIP_DEMUX" = 0 ]]; then
    mkdir -p "$DEMUX_FASTQ_DIR"
  fi
}

prepare_workspace_for_mapping() {
  reset_mapping_outputs
  prepare_workspace_dirs
  mkdir -p "$COUNTS_DIR" "$BAM_DIR"
}

prepare_workspace_for_bigwig() {
  reset_bigwig_outputs
  prepare_workspace_dirs
  mkdir -p "$BIGWIG_DIR" "$UMI_BAM_DIR"
}

finalize_outputs() {
  local count_file

  # 診断ログは既に log/ に直接書かれている。scratch の tmp/ は破棄する
  # (run_primary_command の RETURN trap でも掃除されるが、成功時はここで確実に消す)。
  safe_rm_rf "$TMP_DIR"

  shopt -s nullglob
  for count_file in "${COUNTS_DIR:?}/"*; do
    pigz -f -p 1 "$count_file"
  done
  shopt -u nullglob
}
