#!/usr/bin/env bash

# 役割:
#   pic の shell script で共通して使う補助関数を定義する。
# 入力:
#   シェル環境および関数引数。
# 出力:
#   共通パスを設定し、必要に応じてメッセージを表示する。

declare -g PIC_DEBUG=0

enable_debug_mode() {
  PIC_DEBUG=1
}

debug_enabled() {
  [[ "${PIC_DEBUG:-0}" = 1 ]]
}

debug_log() {
  debug_enabled || return 0
  echo "[DEBUG] $*" >&2
}

debug_kv() {
  local key="$1"
  local value="${2:-}"
  debug_log "${key}=${value}"
}

debug_command() {
  debug_enabled || return 0
  local rendered=""
  local arg
  for arg in "$@"; do
    if [[ -n "$rendered" ]]; then
      rendered+=" "
    fi
    rendered+="$(printf '%q' "$arg")"
  done
  debug_log "CMD=${rendered}"
}

wait_for_available_slot() {
  local max_jobs="$1"

  while (( $(jobs -pr | wc -l) >= max_jobs )); do
    wait -n || return 1
  done
}

run_numbered_jobs_in_parallel() {
  local max_jobs="$1"
  local total_jobs="$2"
  local worker_function="$3"
  local job_id

  for job_id in $(seq "$total_jobs"); do
    wait_for_available_slot "$max_jobs" || return 1
    "$worker_function" "$job_id" &
  done

  wait
}

run_list_jobs_in_parallel() {
  local max_jobs="$1"
  local worker_function="$2"
  shift 2

  local job_arg
  for job_arg in "$@"; do
    wait_for_available_slot "$max_jobs" || return 1
    "$worker_function" "$job_arg" &
  done

  wait
}

pic_resolve_lib_dir() {
  if [[ -z "${PIC_LIB:-}" ]]; then
    PIC_LIB="${HOME}/local/lib/pic"
    export PIC_LIB
    debug_log "PIC_LIB was not set. Using default: ${PIC_LIB}"
  fi
  printf "%s\n" "${PIC_LIB}"
}

require_runtime_lib_dir() {
  local lib_dir
  lib_dir="$(pic_resolve_lib_dir)" || exit 1

  if [[ ! -d "$lib_dir" ]]; then
    mkdir -p "$lib_dir" || {
      handle_error "Failed to create runtime lib directory: ${lib_dir}"
      exit 1
    }
    debug_log "Created runtime lib directory: ${lib_dir}"
  fi
}

resolve_existing_path() {
  local input_path="$1"
  if [[ -d "$input_path" ]]; then
    (cd "$input_path" && pwd)
    return
  fi

  local parent_dir base_name
  parent_dir="$(dirname "$input_path")"
  base_name="$(basename "$input_path")"
  (cd "$parent_dir" && printf "%s/%s\n" "$(pwd)" "$base_name")
}

handle_error() {
  echo "[ERROR] $1" >&2
}

log_info() {
  echo "[INFO] $1"
}
