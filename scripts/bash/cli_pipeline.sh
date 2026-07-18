#!/usr/bin/env bash

# 役割:
#   mapping -> deseq2 -> enrich を 1 コマンドで一括実行する `pic all`。
# 入力:
#   --sample-sheet / --raw-fastq-dir / --run-name / --out-dir の 4 オプションのみ。
# 出力:
#   <out-dir>/            : mapping 成果物 (counts, bam, bw ...) と report_*.html / sample_sheet.tsv
#   <out-dir>/summary/    : 集計/中間テーブル (mapping_sum, deftable, aggregate_profile)
#   <out-dir>/deseq2/      : DESeq2 結果 (ファイル名 suffix <run>_<genome> で genome 識別)
#   <out-dir>/enrich/      : ORA/GSEA 結果 (同上)
# 注記:
#   細かい調整 (--fdr / --methods 等) が必要な場合は deseq2 / enrich を
#   個別に実行する。pic all はシンプル実行専用なので余計なオプションは持たない。

# サンプルシートに含まれる全 genome について生成された deftable を列挙し、
# ファイル名 (deftable_<run_name>_<genome>.tsv) から genome 名を取り出す。
pipeline_list_genomes() {
  local f base genome
  for f in "${DEFTABLE_DIR}"/deftable_"${RUN_NAME}"_*.tsv; do
    [[ -e "$f" ]] || continue
    base="$(basename "$f")"
    base="${base%.tsv}"
    genome="${base#deftable_"${RUN_NAME}"_}"
    printf "%s\n" "$genome"
  done
}

run_pipeline_subcommand() {
  if [[ $# -eq 0 ]]; then
    show_help_for_subcommand all
    return 0
  fi

  local sample_sheet="" raw_fastq_dir="" run_name="" out_dir="" demux_fastq_dir=""
  local run_name_set=0 no_report=0

  while [[ $# -gt 0 ]]; do
    case "$1" in
    --sample-sheet)
      sample_sheet="${2:-}"
      shift 2
      ;;
    --raw-fastq-dir)
      raw_fastq_dir="${2:-}"
      shift 2
      ;;
    --demux-fastq-dir)
      demux_fastq_dir="${2:-}"
      shift 2
      ;;
    --run-name)
      run_name="${2:-}"
      run_name_set=1
      shift 2
      ;;
    --out-dir)
      out_dir="${2:-}"
      shift 2
      ;;
    --no-report)
      no_report=1
      shift
      ;;
    --help)
      show_help_for_subcommand all
      return 0
      ;;
    *)
      handle_error "Unknown option: $1"
      handle_error "pic all のオプション: --sample-sheet / --raw-fastq-dir / --demux-fastq-dir / --run-name / --out-dir / --no-report"
      handle_error "詳細な調整が必要な場合は deseq2 / enrich を個別に実行してください。"
      exit 1
      ;;
    esac
  done

  # mapping 用の共通オプションをパースし、グローバル設定を確定する。
  # subcommand を "mapping" にすることで --run-name / --demux-fastq-dir が許可される。
  local mapping_args=()
  mapping_args+=(--out-dir "$out_dir")
  mapping_args+=(--sample-sheet "$sample_sheet")
  if [[ -n "$demux_fastq_dir" ]]; then
    # 分類済み/demux 済み FASTQ を入力に使う (demux をスキップ)。
    mapping_args+=(--demux-fastq-dir "$demux_fastq_dir")
  else
    mapping_args+=(--raw-fastq-dir "$raw_fastq_dir")
  fi
  if [[ "$run_name_set" = 1 ]]; then
    mapping_args+=(--run-name "$run_name")
  fi

  parse_pic_common_long_options "mapping" "${mapping_args[@]}"

  # --- Step 1: mapping ---
  log_info "pic all: [1/5] mapping を実行します"
  run_primary_command

  # --- genome 列挙 ---
  local genomes=()
  local g
  while IFS= read -r g; do
    [[ -n "$g" ]] && genomes+=("$g")
  done < <(pipeline_list_genomes)

  if [[ "${#genomes[@]}" -eq 0 ]]; then
    handle_error "deftable が見つかりません: ${DEFTABLE_DIR}/deftable_*_${RUN_NAME}.tsv"
    handle_error "mapping が DESeq2 用 deftable を生成したか確認してください。"
    exit 1
  fi

  log_info "pic all: 解析対象 genome: ${genomes[*]}"

  # --- Step 2: aggregate (遺伝子全体 TSS→TES 集計) ---
  log_info "pic all: [2/5] aggregate を実行します"
  if ! run_aggregate_subcommand --out-dir "$OUTPUT_DIR" --run-name "$RUN_NAME" --threads "$THREADS"; then
    handle_error "pic all: aggregate に失敗しました (解析は継続します)"
  fi

  # --- Step 3/4: genome ごとに deseq2 -> enrich ---
  local project deftable deseq2_out stats_csv enrich_out deg_clusters
  for g in "${genomes[@]}"; do
    project="${RUN_NAME}_${g}"
    deftable="${DEFTABLE_DIR}/deftable_${project}.tsv"
    deseq2_out="${OUTPUT_DIR}/deseq2"

    log_info "pic all: [3/5] deseq2 を実行します (genome=${g})"
    Rscript "${PIC_R_DIR}/cmd_run_deseq2.R" \
      --deftable "$deftable" \
      --count-dir "$COUNTS_DIR" \
      --genome "$g" \
      --out-dir "$deseq2_out"

    stats_csv="${deseq2_out}/stats_${project}.csv"
    deg_clusters="${deseq2_out}/DEG/DEGCluster/DEGCluster_gene_for_ora_${project}.csv"
    enrich_out="${OUTPUT_DIR}/enrich"

    log_info "pic all: [4/5] enrich を実行します (genome=${g})"
    if [[ -f "$deg_clusters" ]]; then
      Rscript "${PIC_R_DIR}/cmd_run_enrich.R" \
        --stats "$stats_csv" \
        --genome "$g" \
        --out-dir "$enrich_out" \
        --deg-clusters "$deg_clusters"
    else
      Rscript "${PIC_R_DIR}/cmd_run_enrich.R" \
        --stats "$stats_csv" \
        --genome "$g" \
        --out-dir "$enrich_out"
    fi
  done

  # --- Step 5: HTML レポート ---
  if [[ "$no_report" = 1 ]]; then
    log_info "pic all: [5/5] --no-report によりレポート生成をスキップします"
  else
    log_info "pic all: [5/5] HTML レポートを生成します"
    if ! Rscript "${PIC_R_DIR}/cmd_build_report.R" --out-dir "$OUTPUT_DIR"; then
      handle_error "pic all: レポート生成に失敗しました (解析結果は ${OUTPUT_DIR} に保存済み)"
    fi
  fi

  log_info "pic all: 完了しました (out-dir=${OUTPUT_DIR})"
}
