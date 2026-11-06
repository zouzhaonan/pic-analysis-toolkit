#!/usr/bin/env bash

# 役割:
#   mapping -> aggregate / deseq2 / enrich -> report を 1 コマンドで一括実行する `pic all`。
#   mapping の後は依存関係の順に、待たなくてよい部分を同時に走らせる (pipeline_downstream_genome)。
# 入力:
#   --sample-sheet / --raw-fastq-dir / --run-name / --out-dir が基本。
#   --threads / --jobs は mapping に、--threads は aggregate / deseq2 / enrich / report にも渡す。
#   同時に走る部分はそれぞれ --threads を受ける (メモリの枠は環境変数 PIC_MEM_GB)。
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

# pic all が背後で走らせた R を止める (pic all が止められたとき・失敗したとき)。
# R の fork した子は、親がいなくなると自分で終わる (analysis_runtime.R)。
pipeline_kill_children() {
  local p
  for p in "${PIPELINE_PIDS[@]:-}"; do
    [[ -n "$p" ]] && kill -TERM "$p" 2>/dev/null || true
  done
}

# 1 genome の deseq2 と enrich を依存関係の順に流す。失敗したら 1 を返す。
pipeline_downstream_genome() {
  local g="$1"
  local project="${RUN_NAME}_${g}"
  local deftable="${DEFTABLE_DIR}/deftable_${project}.tsv"
  local deseq2_out="${OUTPUT_DIR}/deseq2"
  local enrich_out="${OUTPUT_DIR}/enrich"
  local stats_csv="${deseq2_out}/stats_${project}.csv"
  # 平坦化後の DEG サブフォルダは suffix 付き (DEG_<proj>/DEGCluster_<proj>/)。
  local deg_clusters="${deseq2_out}/DEG_${project}/DEGCluster_${project}/DEGCluster_gene_for_ora_${project}.csv"
  local de_args=(--deftable "$deftable" --count-dir "$COUNTS_DIR" --genome "$g" --out-dir "$deseq2_out" --threads "$THREADS")
  local fit_pid rlog_pid gsea_pid rc=0

  log_info "pic all: deseq2 (fit と rlog を同時に) を始めます (genome=${g})"
  Rscript "${PIC_R_DIR}/cmd_run_deseq2.R" "${de_args[@]}" --part fit &
  fit_pid=$!
  Rscript "${PIC_R_DIR}/cmd_run_deseq2.R" "${de_args[@]}" --part rlog &
  rlog_pid=$!
  PIPELINE_PIDS+=("$fit_pid" "$rlog_pid")

  if ! wait "$fit_pid"; then
    handle_error "pic all: deseq2 --part fit に失敗しました (genome=${g})"
    kill -TERM "$rlog_pid" 2>/dev/null || true
    wait "$rlog_pid" 2>/dev/null || true
    return 1
  fi

  log_info "pic all: enrich (GSEA) を始めます (genome=${g})"
  Rscript "${PIC_R_DIR}/cmd_run_enrich.R" --stats "$stats_csv" --genome "$g" --out-dir "$enrich_out" \
    --threads "$THREADS" --analyses gsea &
  gsea_pid=$!
  PIPELINE_PIDS+=("$gsea_pid")

  if ! wait "$rlog_pid"; then
    handle_error "pic all: deseq2 --part rlog に失敗しました (genome=${g})"
    rc=1
  elif ! Rscript "${PIC_R_DIR}/cmd_run_deseq2.R" "${de_args[@]}" --part cluster; then
    handle_error "pic all: deseq2 --part cluster に失敗しました (genome=${g})"
    rc=1
  elif [[ -f "$deg_clusters" ]]; then
    log_info "pic all: enrich (ORA) を始めます (genome=${g})"
    if ! Rscript "${PIC_R_DIR}/cmd_run_enrich.R" --stats "$stats_csv" --genome "$g" --out-dir "$enrich_out" \
      --threads "$THREADS" --analyses ora --deg-clusters "$deg_clusters"; then
      handle_error "pic all: enrich (ORA) に失敗しました (genome=${g})"
      rc=1
    fi
  else
    log_info "pic all: DEG クラスタが無いので ORA は行いません (genome=${g})"
  fi

  if [[ "$rc" != 0 ]]; then
    kill -TERM "$gsea_pid" 2>/dev/null || true
    wait "$gsea_pid" 2>/dev/null || true
    return 1
  fi
  if ! wait "$gsea_pid"; then
    handle_error "pic all: enrich (GSEA) に失敗しました (genome=${g})"
    return 1
  fi
  return 0
}

# mapping の後 (aggregate / deseq2 / enrich / report) を流す。pic all と試験の両方から呼ぶ。
#   pipeline_run_downstream <no_aggregate 0|1> <no_report 0|1> <genome...>
# OUTPUT_DIR・RUN_NAME・DEFTABLE_DIR・COUNTS_DIR・THREADS・PIC_R_DIR を使う。
pipeline_run_downstream() {
  # --- Step 2-5: aggregate / deseq2 / enrich / report を依存関係の順に流す ---
  # 1 本で回る段が続かないよう、待たなくてよいものは同時に走らせる:
  #   aggregate ............................................. (ほかと独立。失敗しても続ける)
  #   deseq2 --part fit ──┬─> enrich --analyses gsea
  #   deseq2 --part rlog ─┴─> deseq2 --part cluster ─> enrich --analyses ora
  #   すべて済んだら report
  # どの部分も --threads を受け、worker は仕事の数とメモリの枠 (PIC_MEM_GB) でだけ減らす。
  # 出力は 1 本ずつ順に回したときと同じ。
  local no_aggregate="$1" no_report="$2"
  shift 2
  local genomes=("$@") g
  PIPELINE_PIDS=()
  trap 'pipeline_kill_children' TERM INT HUP

  local agg_pid=""
  if [[ "$no_aggregate" = 1 ]]; then
    log_info "pic all: --no-aggregate により aggregate をスキップします"
  else
    log_info "pic all: aggregate を始めます"
    run_aggregate_subcommand --out-dir "$OUTPUT_DIR" --run-name "$RUN_NAME" --threads "$THREADS" &
    agg_pid=$!
    PIPELINE_PIDS+=("$agg_pid")
  fi

  local failed=0
  for g in "${genomes[@]}"; do
    pipeline_downstream_genome "$g" || failed=1
    [[ "$failed" = 0 ]] || break
  done

  if [[ -n "$agg_pid" ]]; then
    if ! wait "$agg_pid"; then
      handle_error "pic all: aggregate に失敗しました (解析は継続します)"
    fi
  fi
  if [[ "$failed" != 0 ]]; then
    pipeline_kill_children
    handle_error "pic all: deseq2 / enrich に失敗しました (out-dir=${OUTPUT_DIR})"
    exit 1
  fi
  trap - TERM INT HUP

  if [[ "$no_report" = 1 ]]; then
    log_info "pic all: --no-report によりレポート生成をスキップします"
  else
    log_info "pic all: HTML レポートを生成します"
    if ! Rscript "${PIC_R_DIR}/cmd_build_report.R" --out-dir "$OUTPUT_DIR" --threads "$THREADS"; then
      handle_error "pic all: レポート生成に失敗しました (解析結果は ${OUTPUT_DIR} に保存済み)"
    fi
  fi
}

run_pipeline_subcommand() {
  if [[ $# -eq 0 ]]; then
    show_help_for_subcommand all
    return 0
  fi

  local sample_sheet="" raw_fastq_dir="" run_name="" out_dir="" demux_fastq_dir=""
  local run_name_set=0 no_report=0 no_aggregate=0 barcode_first=0
  local hisat2_very_sensitive=0 hisat2_score_min=""
  # 既定は pic mapping と同じ (threads=PIC_DEFAULT_THREADS, jobs=1)。
  # 指定がなければ mapping に渡さず、mapping 側の既定に任せる。
  local threads="" jobs=""

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
    --no-aggregate)
      no_aggregate=1
      shift
      ;;
    --barcode-first)
      barcode_first=1
      shift
      ;;
    --hisat2-very-sensitive)
      hisat2_very_sensitive=1
      shift
      ;;
    --hisat2-score-min)
      hisat2_score_min="${2:-}"
      shift 2
      ;;
    --threads)
      threads="${2:-}"
      shift 2
      ;;
    --jobs)
      jobs="${2:-}"
      shift 2
      ;;
    --help)
      show_help_for_subcommand all
      return 0
      ;;
    *)
      handle_error "Unknown option: $1"
      handle_error "pic all のオプション: --sample-sheet / --raw-fastq-dir / --demux-fastq-dir / --run-name / --out-dir / --threads / --jobs / --no-report / --no-aggregate / --hisat2-score-min / --hisat2-very-sensitive"
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
  if [[ "$hisat2_very_sensitive" = 1 ]]; then
    mapping_args+=(--hisat2-very-sensitive)
  fi
  if [[ -n "$hisat2_score_min" ]]; then
    mapping_args+=(--hisat2-score-min "$hisat2_score_min")
  fi
  if [[ -n "$threads" ]]; then
    mapping_args+=(--threads "$threads")
  fi
  if [[ -n "$jobs" ]]; then
    mapping_args+=(--jobs "$jobs")
  fi
  if [[ "$barcode_first" = 1 ]]; then
    # read1 が barcode(6) + UMI(6) 並びのライブラリ (既定は UMI + barcode)
    mapping_args+=(--barcode-first)
  fi

  # threads / jobs の検証もここで行われ、確定値は THREADS / MAPPING_JOBS に入る。
  # 以降の aggregate / enrich は THREADS を使う。
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

  pipeline_run_downstream "$no_aggregate" "$no_report" "${genomes[@]}"

  log_info "pic all: 完了しました (out-dir=${OUTPUT_DIR})"
}
