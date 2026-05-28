#!/usr/bin/env bash

# 役割:
#   pic のオプションに使う初期値をまとめ、実行時に使う
#   ディレクトリや参照先を設定する。

declare -gr PIC_DEFAULT_THREADS=8
declare -gr PIC_DEFAULT_SKIP_DEMUX=0
declare -gr PIC_DEFAULT_SIMULATION_MODE=0
declare -gr PIC_DEFAULT_HISAT2_VERY_SENSITIVE=0
declare -gr PIC_DEFAULT_OUTPUT_DIR=""
declare -gr PIC_DEFAULT_SAMPLE_SHEET=""
declare -gr PIC_DEFAULT_DEMUX_FASTQ_DIR=""
declare -gr PIC_DEFAULT_RAW_FASTQ_DIR=""

pic_default_run_name() {
  local chars='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
  local result=""
  local i idx

  for i in {1..8}; do
    idx=$((RANDOM % ${#chars}))
    result+="${chars:idx:1}"
  done

  printf "%s\n" "$result"
}

# 役割:
#   pic を動かすときに使うディレクトリや参照先をまとめて設定する。
# 入力:
#   出力ディレクトリ、sample sheet、thread 数、各フラグ、
#   run 名、必要なら入力に使う demux 済み FASTQ ディレクトリ、
#   hisat2 の追加フラグ。
# 出力:
#   グローバル変数を設定する。

init_config() {
  local output_dir="$1"
  local sample_sheet="$2"
  local threads="$3"
  local skip_demux="$4"
  local simulation_mode="$5"
  local run_name="$6"
  local requested_demux_fastq_dir="${7:-}"
  local hisat2_very_sensitive="$8"
  local raw_fastq_dir="${9:-}"
  local demux_layout="${10:-nested}"

  declare -gr SKIP_DEMUX="$skip_demux"
  declare -gr SIMULATION_MODE="$simulation_mode"
  declare -gr HISAT2_VERY_SENSITIVE="$hisat2_very_sensitive"
  declare -gr LIB_DIR="$(pic_resolve_lib_dir)"
  declare -gr PIC_GENOME_MAP_FILE="${LIB_DIR}/register/genome_map.tsv"

  declare -gr PIC_DIR="$PIC_ROOT"
  declare -gr OUTPUT_DIR="$output_dir"
  declare -gr RUN_NAME="$run_name"
  declare -gr BAM_DIR="$output_dir/bam"
  declare -gr UMI_BAM_DIR="$output_dir/umi_bam"
  declare -gr BIGWIG_DIR="$output_dir/bw"
  declare -gr TMP_DIR="$output_dir/tmp"
  declare -gr COUNTS_DIR="$output_dir/counts"
  declare -gr LOG_DIR="$output_dir/log"

  declare -gr SAMPLE_SHEET="$sample_sheet"
  declare -gr RAW_FASTQ_DIR="$raw_fastq_dir"
  declare -gr MAP_SUMMARY="${OUTPUT_DIR}/mapping_sum__${run_name}.tsv"
  declare -gr DEFTABLE_DIR="${OUTPUT_DIR}"
  declare -gr META_FOR_CHUNK="${TMP_DIR}/meta_for_chunk"
  declare -gr FORMATTED_SAMPLE_SHEET="${TMP_DIR}/formatted_sample"
  declare -gr JOB_QUEUE="${TMP_DIR}/job_queue"

  declare -gr FORMAT_SAMPLE_AWK="${PIC_DIR}/scripts/awk/format_sample.awk"
  declare -gr SPLIT_FASTQ_AWK="${PIC_DIR}/scripts/awk/split_to_chunks.awk"
  declare -gr MAP_SUMMARY_AWK="${PIC_DIR}/scripts/awk/mapping_sum.awk"
  declare -gr PREPARE_JOB_QUEUE_AWK="${PIC_DIR}/scripts/awk/map_prepare_job_queue.awk"
  declare -gr GENERATE_DEFTABLES_AWK="${PIC_DIR}/scripts/awk/summarize_generate_deftables.awk"
  declare -gr EXTRACT_SAMPLE_SHEET_GENOMES_AWK="${PIC_DIR}/scripts/awk/validate_extract_sample_sheet_genomes.awk"
  declare -gr PLOT_SIM_SCRIPT="${PIC_DIR}/scripts/r/pic_plot_sim.R"

  local resolved_demux_fastq_dir="${output_dir}/demux"
  if [[ "$demux_layout" = "flat" ]]; then
    resolved_demux_fastq_dir="$output_dir"
  fi
  if [[ "$SKIP_DEMUX" = 1 ]]; then
    # このモードでは、既に demux 済みの FASTQ ディレクトリを入力に使う。
    resolved_demux_fastq_dir="$requested_demux_fastq_dir"
  fi

  declare -gr DEMUX_FASTQ_DIR="$resolved_demux_fastq_dir"
  declare -gr TOTAL_READS="${DEMUX_FASTQ_DIR}/total_reads"
  declare -gr THREADS="$threads"
}

init_prepare_config() {
  local genome_name="$1"
  local fasta_path="$2"
  local gtf_path="$3"
  local fasta_source="$4"
  local gtf_source="$5"
  local threads="$6"

  declare -gr LIB_DIR="$(pic_resolve_lib_dir)"
  declare -gr PIC_GENOME_MAP_FILE="${LIB_DIR}/register/genome_map.tsv"
  declare -gr PREPARE_GENOME="$genome_name"
  declare -gr PREPARE_FASTA="$fasta_path"
  declare -gr PREPARE_GTF="$gtf_path"
  declare -gr PREPARE_FASTA_SOURCE="$fasta_source"
  declare -gr PREPARE_GTF_SOURCE="$gtf_source"
  declare -gr PREPARE_THREADS="$threads"
}
