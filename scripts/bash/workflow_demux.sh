#!/usr/bin/env bash

# demux workflow:
#   FASTQ の組 (<prefix>_R1_001 / _R2_001) ごとに 1 回だけ読み、barcode でサンプルへ振り分ける
#   (demux_single_prefix)。prefix ごとの出力 (part) を、最後にサンプルごとにつなぐ
#   (collect_demux_parts)。
#
#   part:    <dir>/<sample>.part.<prefix 名>.fastq.gz   (読みが 1 つも無いサンプルは作らない)
#   済んだ印: <dir>/total_reads.part.<prefix 名>         (sample<TAB>読みの数。最後に書く)
#   まとめ:  <demux>/<sample>.fastq.gz と <demux>/total_reads (今までの pic demux と同じ)
#
#   <prefix 名> は sample sheet の 1 列目の basename。
#   旧版 (seqkit range で THREADS 個に分けて読み直し、あとで cat | pigz) と出力の中身は同じ。
#   数え直し (gzcat | wc -l)・読み直し・まとめの待ち合わせが無くなる。

# sample sheet の本体の行 (ヘッダーと CR を除き、1 列目が空の行を飛ばす)。
demux_sheet_rows() {
  tr -d "\015" <"$SAMPLE_SHEET" | awk -F"\t" -v OFS="\t" '
    NR == 1 && tolower($1) == "fastq_prefix" && tolower($2) == "genome" && tolower($3) == "barcode" && tolower($4) == "sample" && tolower($5) == "group" { next }
    $1 != "" { print }'
}

demux_prefix_name() {
  basename "$1"
}

# sample sheet の 1 列目の値 (sheet の順、重複なし)。
list_demux_prefixes() {
  demux_sheet_rows | awk -F"\t" '!seen[$1]++ { print $1 }'
}

demux_raw_fastq_prefix() {
  local prefix="$1"
  if [[ "$prefix" == /* ]]; then
    printf "%s\n" "$prefix"
  else
    printf "%s/%s\n" "$RAW_FASTQ_DIR" "$prefix"
  fi
}

demux_part_marker() {
  local part_dir="$1"
  local name="$2"
  printf "%s/total_reads.part.%s\n" "$part_dir" "$name"
}

# 1 列目の値が違うのに basename が同じ prefix があると、part の名前がぶつかる。
validate_demux_prefix_names() {
  local dup
  dup="$(list_demux_prefixes | while IFS= read -r p; do basename "$p"; done | sort | uniq -d)"
  if [[ -n "$dup" ]]; then
    handle_error "Different fastq_prefix values share the same name: ${dup//$'\n'/, }"
    exit 1
  fi
}

# 1 つの prefix を振り分け、part と済んだ印を part_dir に書く。
# 引数の prefix は sheet の 1 列目の値か、その basename。threads は pic-demux に渡すスレッド数。
# 戻り値: 0 = 済んだ (もともと済んでいた場合も含む)、1 = 失敗。
demux_single_prefix() {
  local prefix_arg="$1"
  local part_dir="$2"
  local threads="${3:-$THREADS}"
  local rows sheet_prefix name fq_prefix barcodes samples marker suffix
  local s rc=0
  local barcode_first_opt=()

  rows="$(demux_sheet_rows | awk -F"\t" -v OFS="\t" -v p="$prefix_arg" '{ n = $1; sub(/.*\//, "", n) } $1 == p || n == p')"
  if [[ -z "$rows" ]]; then
    handle_error "fastq_prefix not found in sample sheet: ${prefix_arg}"
    return 1
  fi
  sheet_prefix="$(printf "%s\n" "$rows" | head -n 1 | cut -f1)"
  name="$(demux_prefix_name "$sheet_prefix")"
  fq_prefix="$(demux_raw_fastq_prefix "$sheet_prefix")"
  marker="$(demux_part_marker "$part_dir" "$name")"
  suffix=".part.${name}.fastq.gz"

  if [[ -f "$marker" ]]; then
    log_info "Skipping ${name}: already done (${marker})"
    return 0
  fi

  barcodes="$(printf "%s\n" "$rows" | awk -F"\t" '{ printf "%s%s", (NR > 1 ? "," : ""), $3 }')"
  samples="$(printf "%s\n" "$rows" | awk -F"\t" '{ printf "%s%s", (NR > 1 ? "," : ""), $4 }')"

  if ! command -v "${PIC_DEMUX_BIN:-pic-demux}" >/dev/null 2>&1; then
    handle_error "pic-demux not found. Build it with: bash ${PIC_ROOT:-<toolkit>}/src/pic-demux/install.sh"
    return 1
  fi

  local fq
  for fq in "${fq_prefix}_R1_001.fastq.gz" "${fq_prefix}_R2_001.fastq.gz"; do
    if [[ ! -f "$fq" ]]; then
      handle_error "FASTQ not found: ${fq}"
      return 1
    fi
  done

  # 前に途中で止まった分を消す (このサンプルとこの prefix の組だけ)。
  while IFS= read -r s; do
    rm -f "${part_dir}/${s}${suffix}" "${part_dir}/${s}${suffix}.tmp"
  done < <(printf "%s\n" "$rows" | cut -f4 | sort -u)
  rm -f "${marker}.tmp" "${marker}.nread.tmp"
  mkdir -p "$part_dir"

  log_info "Demultiplexing ${name}"
  debug_log "STEP=demux_single_prefix prefix=${name}"
  debug_kv "input.fastq_prefix" "$fq_prefix"
  debug_kv "arg.barcodes" "$barcodes"
  debug_kv "arg.samples" "$samples"

  # pic-demux (src/pic-demux): R1 と R2 を別々のスレッドで解凍して突き合わせ、サンプルごとの
  # gzip をスレッドで並べて圧縮する。中身は旧 demux_prefix.awk (paste | awk) と同じ。
  # 入力の不正 (gzip の破損・途中切れ、R1 と R2 の行数違い、4 行単位でない) は 0 以外で終わる。
  [[ "${PIC_BARCODE_FIRST:-0}" = 1 ]] && barcode_first_opt=(--barcode-first)
  "${PIC_DEMUX_BIN:-pic-demux}" \
    --r1 "${fq_prefix}_R1_001.fastq.gz" --r2 "${fq_prefix}_R2_001.fastq.gz" \
    --barcodes "$barcodes" --samples "$samples" \
    --out-dir "$part_dir" --suffix "$suffix" --nread-file "${marker}.nread.tmp" \
    --threads "$threads" --level "${PIC_DEMUX_LEVEL:-6}" \
    "${barcode_first_opt[@]}" || rc=$?
  if [[ "$rc" != 0 ]]; then
    handle_error "Demux failed for ${name} (pic-demux status: ${rc})."
    while IFS= read -r s; do
      rm -f "${part_dir}/${s}${suffix}.tmp"
    done < <(printf "%s\n" "$rows" | cut -f4 | sort -u)
    rm -f "${marker}.nread.tmp"
    return 1
  fi

  while IFS= read -r s; do
    if [[ -f "${part_dir}/${s}${suffix}.tmp" ]]; then
      mv "${part_dir}/${s}${suffix}.tmp" "${part_dir}/${s}${suffix}"
    fi
  done < <(printf "%s\n" "$rows" | cut -f4 | sort -u)
  mv "${marker}.nread.tmp" "$marker"
  debug_log "RESULT=demux_single_prefix prefix=${name} status=ok"
}

# 済んでいない prefix 名を出す。
list_demux_prefixes_without_marker() {
  local part_dir="$1"
  local p name
  while IFS= read -r p; do
    name="$(demux_prefix_name "$p")"
    [[ -f "$(demux_part_marker "$part_dir" "$name")" ]] || printf "%s\n" "$name"
  done < <(list_demux_prefixes)
}

# part をサンプルごとにつなぎ (gzip は cat でつなげる)、total_reads を作る。
# 1 つのサンプルが複数の prefix にあるときは、prefix 名の順 (C ロケール) につなぐ。
# 旧版の並び (tmp/chunks/<prefix 名>/ の glob 順) と、普通の prefix 名なら同じ。
collect_demux_parts() {
  local part_dir="$1"
  local dest_dir="$2"
  local missing totals_dir s name dest parts=() names=()

  missing="$(list_demux_prefixes_without_marker "$part_dir")"
  if [[ -n "$missing" ]]; then
    handle_error "Demux not finished yet for (no ${part_dir}/total_reads.part.<prefix>):"
    printf "%s\n" "$missing" >&2
    return 1
  fi

  mkdir -p "$dest_dir"
  while IFS= read -r name; do
    names+=("$name")
  done < <(list_demux_prefixes | while IFS= read -r p; do basename "$p"; done | LC_ALL=C sort -u)

  totals_dir="${TMP_DIR}/demux_totals"
  safe_rm_rf "$totals_dir"
  mkdir -p "$totals_dir"

  while IFS= read -r s; do
    parts=()
    for name in "${names[@]}"; do
      [[ -f "${part_dir}/${s}.part.${name}.fastq.gz" ]] && parts+=("${part_dir}/${s}.part.${name}.fastq.gz")
    done
    awk -F"\t" -v s="$s" '$1 == s { sum += $2 } END { print s "\t" sum + 0 }' \
      "${names[@]/#/${part_dir}/total_reads.part.}" >"${totals_dir}/total_reads.${s}"
    if [[ "${#parts[@]}" -eq 0 ]]; then
      continue
    fi
    dest="$(demux_fastq_path "$s")"
    # part が 1 つだけなら、写さずにハードリンクする (中身は同じ。part は最後に消すので、途中で
    # 止まってもやり直せる)。別のファイルシステムなどでリンクできなければ cat で写す。
    rm -f "${dest}.tmp"
    if [[ "${#parts[@]}" -ne 1 ]] || ! ln "${parts[0]}" "${dest}.tmp" 2>/dev/null; then
      cat "${parts[@]}" >"${dest}.tmp"
    fi
    mv "${dest}.tmp" "$dest"
  done < <(demux_sheet_rows | cut -f4 | sort -u)

  # total_reads の行の並びは旧版と同じ (total_reads.<sample> の glob 順)。
  cat "${totals_dir}"/total_reads.* >"${TOTAL_READS}.tmp"
  mv "${TOTAL_READS}.tmp" "$TOTAL_READS"
  safe_rm_rf "$totals_dir"

  # まとめ終わった part と印を消す。
  for name in "${names[@]}"; do
    rm -f "${part_dir}"/*".part.${name}.fastq.gz" "$(demux_part_marker "$part_dir" "$name")"
  done
  debug_kv "result.total_reads_rows" "$(wc -l <"$TOTAL_READS")"
}

# 全 prefix を並べて振り分け、まとめる (pic demux / pic mapping / pic all の既定)。
# pic-demux 1 本は --threads のぶんのコアを使う (解凍と突き合わせに約 2、残りで圧縮を並べる)。
# 1 本が速くなるのは約 12 コアまで (R2 の解凍と突き合わせが 1 本ずつなので)。1 本に 6 コア以上を
# 割り当て、同時に回す prefix の数は floor(THREADS / 6) (1 以上、prefix の数まで)。回る順番の
# 段の数を変えない範囲で同時の数を減らして揃え、1 本の --threads は THREADS / 同時の数。
# 例: THREADS=40・6 prefix → 6 本 × 6、THREADS=8 → 1 本 × 8 を 6 回、THREADS=24 → 3 本 × 8 を 2 回。
split_fq_by_barcode() {
  local part_dir="${TMP_DIR}/demux_parts"
  local status_dir="${TMP_DIR}/demux_status"
  local prefixes=() p i=0 failed=0 max_jobs job_threads

  validate_demux_prefix_names
  debug_log "STEP=split_fq_by_barcode"
  mkdir -p "$part_dir" "$status_dir"

  while IFS= read -r p; do
    prefixes+=("$p")
  done < <(list_demux_prefixes)

  local n_prefix=${#prefixes[@]} waves
  max_jobs=$(( THREADS / 6 ))
  (( max_jobs > n_prefix )) && max_jobs=$n_prefix
  (( max_jobs < 1 )) && max_jobs=1
  waves=$(( (n_prefix + max_jobs - 1) / max_jobs ))
  (( waves < 1 )) && waves=1
  max_jobs=$(( (n_prefix + waves - 1) / waves ))
  (( max_jobs < 1 )) && max_jobs=1
  job_threads=$(( THREADS / max_jobs ))
  debug_kv "arg.parallel_prefixes" "$max_jobs"
  debug_kv "arg.threads_per_prefix" "$job_threads"

  for p in "${prefixes[@]}"; do
    i=$((i + 1))
    while (( $(jobs -pr | wc -l) >= max_jobs )); do
      wait -n || true
    done
    ( rc=0; demux_single_prefix "$p" "$part_dir" "$job_threads" || rc=$?; echo "$rc" >"${status_dir}/${i}" ) &
  done
  wait || true

  for i in $(seq "${#prefixes[@]}"); do
    if [[ "$(cat "${status_dir}/${i}" 2>/dev/null)" != 0 ]]; then
      failed=1
    fi
  done
  safe_rm_rf "$status_dir"
  if [[ "$failed" = 1 ]]; then
    handle_error "Demux failed for one or more FASTQ prefixes."
    return 1
  fi

  collect_demux_parts "$part_dir" "$DEMUX_FASTQ_DIR" || return 1
  safe_rm_rf "$part_dir"
  debug_log "RESULT=split_fq_by_barcode status=ok"
}

# pic demux (既定 / --prefix / --collect)
run_demux_command() {
  local mode="${PIC_DEMUX_MODE:-all}"

  log_info "Validating inputs"
  validate_inputs

  case "$mode" in
  prefix)
    # 同じ out-dir で別々の prefix を同時に回せるよう、何も消さずに part だけを書く。
    validate_demux_prefix_names
    mkdir -p "$DEMUX_FASTQ_DIR"
    demux_single_prefix "$PIC_DEMUX_PREFIX" "$DEMUX_FASTQ_DIR" || exit 1
    ;;
  collect)
    validate_demux_prefix_names
    mkdir -p "$DEMUX_FASTQ_DIR"
    log_info "Collecting demux parts"
    collect_demux_parts "$DEMUX_FASTQ_DIR" "$DEMUX_FASTQ_DIR" || exit 1
    safe_rm_rf "$TMP_DIR"
    ;;
  *)
    log_info "Preparing demux workspace"
    prepare_workspace_for_demux
    log_info "Running demux by barcode"
    split_fq_by_barcode || exit 1
    safe_rm_rf "$TMP_DIR"
    ;;
  esac
  debug_log "RESULT=run_demux_command mode=${mode} status=ok"
}
