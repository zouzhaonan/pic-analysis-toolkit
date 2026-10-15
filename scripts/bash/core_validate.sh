#!/usr/bin/env bash

# 役割:
#   pic を動かす前に、設定値と入力ファイル・入力フォルダを確認する。
# 入力:
#   init_config で設定されたグローバル変数。
# 出力:
#   なし。問題があればエラーメッセージを出して終了する。

validate_common_arguments() {
  local threads="$1"
  local output_dir="$2"
  local sample_sheet="$3"
  local skip_demux="$4"
  local demux_fastq_dir="$5"
  local run_name="$6"
  local raw_fastq_dir="$7"

  if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
    handle_error "Threads must be a positive integer: ${threads}"
    exit 1
  fi

  if [[ -z "$output_dir" || -z "$sample_sheet" ]]; then
    handle_error "Both output directory and sample sheet are required."
    exit 1
  fi

  if [[ "$skip_demux" = 1 && -z "$demux_fastq_dir" ]]; then
    handle_error "Demux FASTQ directory is required when --demux-fastq-dir is used."
    exit 1
  fi

  if [[ "$skip_demux" = 0 && -z "$raw_fastq_dir" ]]; then
    handle_error "Raw FASTQ directory is required unless --demux-fastq-dir is used."
    handle_error "Specify --raw-fastq-dir <dir>."
    exit 1
  fi

  if [[ -z "$run_name" ]]; then
    handle_error "Run name must not be empty."
    exit 1
  fi
}

validate_inputs() {
  if [[ -z "${OUTPUT_DIR:-}" || -z "${SAMPLE_SHEET:-}" ]]; then
    handle_error "Configuration is not initialized."
    exit 1
  fi

  if [[ ! -f "$SAMPLE_SHEET" ]]; then
    handle_error "Sample sheet not found: ${SAMPLE_SHEET}"
    exit 1
  fi

  if [[ ! -s "$SAMPLE_SHEET" ]]; then
    handle_error "Sample sheet is empty: ${SAMPLE_SHEET}"
    exit 1
  fi

  if [[ "$SKIP_DEMUX" = 1 && ! -d "$DEMUX_FASTQ_DIR" ]]; then
    handle_error "Input FASTQ directory not found: ${DEMUX_FASTQ_DIR}"
    exit 1
  fi

  if [[ "$SKIP_DEMUX" = 1 && ! -f "$TOTAL_READS" ]]; then
    handle_error "total_reads not found in input FASTQ directory: ${TOTAL_READS}"
    exit 1
  fi

  if [[ "$SKIP_DEMUX" = 0 && ! -d "$RAW_FASTQ_DIR" ]]; then
    handle_error "FASTQ prefix directory not found: ${RAW_FASTQ_DIR}"
    exit 1
  fi

  validate_sample_sheet_rows
}

# 同じ sample 名の複数行は「複数 barcode を 1 サンプルに合算」の指定として許すが、
# genome / group が食い違うものと、同じ (fastq_prefix, barcode) の重複は誤記として止める。
# (fastq_prefix, barcode) の重複は raw FASTQ から demux するときだけの問題。
# --demux-fastq-dir では同じ読みを別サンプル (例: barcode 別と合算) に使ってよい。
validate_sample_sheet_rows() {
  local problems
  problems="$(tr -d '\015' <"$SAMPLE_SHEET" | awk -F'\t' -v skip_demux="${SKIP_DEMUX:-0}" '
    NR == 1 && tolower($1) == "fastq_prefix" { next }
    NF < 5 { next }
    {
      key = $1 SUBSEP $3
      if (skip_demux != 1 && (key in pair)) printf "duplicate fastq_prefix/barcode %s/%s (samples %s, %s)\n", $1, $3, pair[key], $4
      pair[key] = $4
      if (($4 in genome) && genome[$4] != $2) printf "sample %s has different genomes (%s, %s)\n", $4, genome[$4], $2
      if (($4 in group) && group[$4] != $5) printf "sample %s has different groups (%s, %s)\n", $4, group[$4], $5
      genome[$4] = $2; group[$4] = $5
    }')"
  if [[ -n "$problems" ]]; then
    handle_error "Invalid sample sheet rows:"$'\n'"${problems}"
    exit 1
  fi
}

validate_sample_sheet_genome_registration() {
  local genomes genome_name
  genomes="$(
    awk -F '\t' -f "$EXTRACT_SAMPLE_SHEET_GENOMES_AWK" "$SAMPLE_SHEET" | sort -u
  )"

  while IFS= read -r genome_name; do
    [[ -z "$genome_name" ]] && continue
    if ! pic_validate_registered_genome "$genome_name"; then
      exit 1
    fi
  done <<< "$genomes"
}

validate_bigwig_inputs() {
  local missing_genomes

  if [[ ! -f "$JOB_QUEUE" ]]; then
    handle_error "job queue not found: ${JOB_QUEUE}"
    handle_error "Run mapping before generating BigWig files."
    exit 1
  fi

  missing_genomes="$(list_missing_chrom_sizes_genomes)"
  if [[ -n "$missing_genomes" ]]; then
    handle_error "chrom size file not found for:"
    printf "%s\n" "$missing_genomes" >&2
    exit 1
  fi
}
