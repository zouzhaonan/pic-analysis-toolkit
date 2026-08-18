#!/usr/bin/env bash

parse_pic_common_long_options() {
  local subcommand="${1:-}"
  shift || true

  local allow_run_name=0
  local allow_demux_fastq=0
  local allow_simulation=0
  local allow_hisat2_very_sensitive=0

  case "$subcommand" in
  mapping)
    allow_run_name=1
    allow_demux_fastq=1
    allow_simulation=1
    allow_hisat2_very_sensitive=1
    ;;
  demux) ;;
  esac

  PIC_WANT_HELP=0
  PIC_OUTPUT_DIR=""
  PIC_SAMPLE_SHEET=""
  PIC_THREADS="$PIC_DEFAULT_THREADS"
  PIC_SKIP_DEMUX="$PIC_DEFAULT_SKIP_DEMUX"
  PIC_SIMULATION_MODE="$PIC_DEFAULT_SIMULATION_MODE"
  PIC_HISAT2_VERY_SENSITIVE="$PIC_DEFAULT_HISAT2_VERY_SENSITIVE"
  PIC_HISAT2_SCORE_MIN="$PIC_DEFAULT_HISAT2_SCORE_MIN"
  PIC_FILTER_PRIMER_READS="$PIC_DEFAULT_FILTER_PRIMER_READS"
  PIC_DEMUX_FASTQ_DIR="$PIC_DEFAULT_DEMUX_FASTQ_DIR"
  PIC_RAW_FASTQ_DIR="$PIC_DEFAULT_RAW_FASTQ_DIR"
  PIC_DEMUX_LAYOUT=""
  PIC_RUN_NAME="$(pic_default_run_name)"

  case "$subcommand" in
  mapping) PIC_DEMUX_LAYOUT="nested" ;;
  demux) PIC_DEMUX_LAYOUT="flat" ;;
  esac

  while [[ $# -gt 0 ]]; do
    case "$1" in
    --out-dir)
      PIC_OUTPUT_DIR="${2:-}"
      shift 2
      ;;
    --sample-sheet)
      PIC_SAMPLE_SHEET="${2:-}"
      shift 2
      ;;
    --threads)
      PIC_THREADS="${2:-}"
      shift 2
      ;;
    --demux-fastq-dir)
      if [[ "$allow_demux_fastq" != 1 ]]; then
        handle_error "--demux-fastq-dir is not supported for ${subcommand}."
        exit 1
      fi
      PIC_DEMUX_FASTQ_DIR="${2:-}"
      PIC_SKIP_DEMUX=1
      shift 2
      ;;
    --raw-fastq-dir)
      PIC_RAW_FASTQ_DIR="${2:-}"
      shift 2
      ;;
    --run-name)
      if [[ "$allow_run_name" != 1 ]]; then
        handle_error "--run-name is not supported for ${subcommand}."
        exit 1
      fi
      PIC_RUN_NAME="${2:-}"
      shift 2
      ;;
    --simulation)
      if [[ "$allow_simulation" != 1 ]]; then
        handle_error "--simulation is not supported for ${subcommand}."
        exit 1
      fi
      PIC_SIMULATION_MODE=1
      shift
      ;;
    --hisat2-very-sensitive)
      if [[ "$allow_hisat2_very_sensitive" != 1 ]]; then
        handle_error "--hisat2-very-sensitive is not supported for ${subcommand}."
        exit 1
      fi
      PIC_HISAT2_VERY_SENSITIVE=1
      shift
      ;;
    --hisat2-score-min)
      if [[ "$allow_hisat2_very_sensitive" != 1 ]]; then
        handle_error "--hisat2-score-min is not supported for ${subcommand}."
        exit 1
      fi
      PIC_HISAT2_SCORE_MIN="${2:-}"
      shift 2
      ;;
    --help)
      PIC_WANT_HELP=1
      return 0
      ;;
    *)
      handle_error "Unknown option: $1"
      exit 1
      ;;
    esac
  done

  validate_common_arguments \
    "$PIC_THREADS" \
    "$PIC_OUTPUT_DIR" \
    "$PIC_SAMPLE_SHEET" \
    "$PIC_SKIP_DEMUX" \
    "$PIC_DEMUX_FASTQ_DIR" \
    "$PIC_RUN_NAME" \
    "$PIC_RAW_FASTQ_DIR"

  init_config \
    "$PIC_OUTPUT_DIR" \
    "$PIC_SAMPLE_SHEET" \
    "$PIC_THREADS" \
    "$PIC_SKIP_DEMUX" \
    "$PIC_SIMULATION_MODE" \
    "$PIC_RUN_NAME" \
    "$PIC_DEMUX_FASTQ_DIR" \
    "$PIC_HISAT2_VERY_SENSITIVE" \
    "$PIC_RAW_FASTQ_DIR" \
    "$PIC_DEMUX_LAYOUT" \
    "$PIC_HISAT2_SCORE_MIN" \
    "$PIC_FILTER_PRIMER_READS"
}

run_primary_family_subcommand() {
  local subcommand="$1"
  shift || true

  if [[ $# -eq 0 ]]; then
    show_help_for_subcommand "$subcommand"
    return 0
  fi

  parse_pic_common_long_options "$subcommand" "$@"
  if [[ "$PIC_WANT_HELP" = 1 ]]; then
    show_help_for_subcommand "$subcommand"
    return 0
  fi

  case "$subcommand" in
  mapping) run_primary_command ;;
  demux) run_demux_command ;;
  *)
    handle_error "Unsupported subcommand: ${subcommand}"
    exit 1
    ;;
  esac
}

run_build_genome_subcommand() {
  if [[ $# -eq 0 ]]; then
    show_help_for_subcommand build-genome
    return 0
  fi

  local genome_name=""
  local fasta_path=""
  local gtf_path=""
  local fasta_source="unspecified"
  local gtf_source="unspecified"
  local threads="$PIC_DEFAULT_THREADS"

  while [[ $# -gt 0 ]]; do
    case "$1" in
    --genome)
      genome_name="${2:-}"
      shift 2
      ;;
    --fasta)
      fasta_path="${2:-}"
      shift 2
      ;;
    --gtf)
      gtf_path="${2:-}"
      shift 2
      ;;
    --fasta-source)
      fasta_source="${2:-}"
      shift 2
      ;;
    --gtf-source)
      gtf_source="${2:-}"
      shift 2
      ;;
    --threads)
      threads="${2:-}"
      shift 2
      ;;
    --help)
      show_help_for_subcommand build-genome
      return 0
      ;;
    *)
      handle_error "Unknown option: $1"
      exit 1
      ;;
    esac
  done

  if [[ -z "$genome_name" || -z "$fasta_path" || -z "$gtf_path" ]]; then
    handle_error "build-genome requires --genome, --fasta and --gtf."
    exit 1
  fi
  if [[ ! -f "$fasta_path" ]]; then
    handle_error "FASTA not found: ${fasta_path}"
    exit 1
  fi
  if [[ ! -f "$gtf_path" ]]; then
    handle_error "GTF not found: ${gtf_path}"
    exit 1
  fi

  init_prepare_config \
    "$genome_name" \
    "$(resolve_existing_path "$fasta_path")" \
    "$(resolve_existing_path "$gtf_path")" \
    "$fasta_source" \
    "$gtf_source" \
    "$threads"

  debug_log "STEP=run_build_genome_subcommand"
  debug_kv "arg.genome" "$genome_name"
  debug_kv "input.fasta" "$fasta_path"
  debug_kv "input.gtf" "$gtf_path"
  debug_kv "arg.fasta_source" "$fasta_source"
  debug_kv "arg.gtf_source" "$gtf_source"
  debug_kv "arg.threads" "$threads"

  run_build_genome_command
}
