#!/usr/bin/env bash

# mapping workflow:
# sample sheet formatting -> optional demux -> mapping/count -> UMI count -> summary -> bigwig

format_sample_sheet() {
  debug_log "STEP=format_sample_sheet"
  debug_kv "input.sample_sheet" "$SAMPLE_SHEET"
  debug_kv "input.raw_fastq_dir" "${RAW_FASTQ_DIR:-}"
  debug_kv "input.awk" "$FORMAT_SAMPLE_AWK"
  debug_kv "output.formatted_sample_sheet" "$FORMATTED_SAMPLE_SHEET"
  debug_kv "output.meta_for_chunk" "$META_FOR_CHUNK"
  {
    printf "%s\t%s\t%s\t%s\t%s\n" "fastq_prefix" "genome" "barcode" "sample" "group"
    tr -d "\015" <"$SAMPLE_SHEET" | awk -F"\t" 'NR==1{
      if (tolower($1)=="fastq_prefix" && tolower($2)=="genome" && tolower($3)=="barcode" && tolower($4)=="sample" && tolower($5)=="group") next
    } {print}'
  } >"$FORMATTED_SAMPLE_SHEET"

  if [[ "$SKIP_DEMUX" = 1 ]]; then
    : >"$META_FOR_CHUNK"
  else
    awk -F"\t" -v OFS="\t" -v threads="$THREADS" -v raw_fastq_dir="$RAW_FASTQ_DIR" \
      -f "$FORMAT_SAMPLE_AWK" "$FORMATTED_SAMPLE_SHEET" >"$META_FOR_CHUNK"
  fi
  debug_kv "result.meta_rows" "$(wc -l <"$META_FOR_CHUNK")"
}

chunk_single_fq() {
  local job_id="$1"
  local chunk fq_prefix range barcodes samples
  local fq_files chunk_tmp_dir

  read -r chunk fq_prefix range barcodes samples < <(
    sed -n "${job_id}p" "$META_FOR_CHUNK"
  )

  fq_files=("${fq_prefix}_R1_001.fastq.gz" "${fq_prefix}_R2_001.fastq.gz")
  chunk_tmp_dir="${TMP_DIR}/chunks/$(basename "$fq_prefix")"
  mkdir -p "$chunk_tmp_dir"

  awk -v barcodes="$barcodes" -v samples="$samples" -v chunk="$chunk" \
    -v tmp_dir="$chunk_tmp_dir" -f "$SPLIT_FASTQ_AWK" <(
      paste \
        <(seqkit range -r "$range" "${fq_files[0]}") \
        <(seqkit range -r "$range" "${fq_files[1]}")
    )
}

chunk_all_input_fq() {
  local total_jobs

  total_jobs="$(wc -l <"$META_FOR_CHUNK")"
  debug_log "STEP=chunk_all_input_fq"
  debug_kv "input.meta_for_chunk" "$META_FOR_CHUNK"
  debug_kv "arg.threads" "$THREADS"
  debug_kv "arg.total_jobs" "$total_jobs"
  run_numbered_jobs_in_parallel "$THREADS" "$total_jobs" chunk_single_fq
  debug_log "RESULT=chunk_all_input_fq done"
}

gather_single_sample_fq() {
  local sample="$1"
  local sample_total_reads_file fastq_chunks nread_chunks

  sample_total_reads_file="${TMP_DIR}/total_reads.${sample}"
  fastq_chunks=("${TMP_DIR}"/chunks/*/"${sample}.fastq".*)
  nread_chunks=("${TMP_DIR}"/chunks/*/"${sample}.nread".*)

  if [[ ! -e "${fastq_chunks[0]}" ]]; then
    printf "%s\t0\n" "$sample" >"$sample_total_reads_file"
    return 0
  fi

  cat "${fastq_chunks[@]}" | pigz -p 1 >"$(demux_fastq_path "$sample")"
  awk -v sample="$sample" '{sum += $1} END {print sample "\t" sum}' \
    "${nread_chunks[@]}" >"$sample_total_reads_file"

  rm -f "${fastq_chunks[@]}" "${nread_chunks[@]}"
}

gather_all_samples_fq() {
  local samples=()
  local sample

  while IFS= read -r sample; do
    samples+=("$sample")
  done < <(tail -n +2 "$FORMATTED_SAMPLE_SHEET" | cut -f4 | sort -u)

  run_list_jobs_in_parallel "$THREADS" gather_single_sample_fq "${samples[@]}"
}

merge_total_reads() {
  cat "${TMP_DIR}"/total_reads.* >"$TOTAL_READS"
  rm -f "${TMP_DIR}"/total_reads.*
}

gather_chunked_fq() {
  debug_log "STEP=gather_chunked_fq"
  debug_kv "input.chunk_dir" "${TMP_DIR}/chunks"
  debug_kv "output.total_reads" "$TOTAL_READS"
  : >"$TOTAL_READS"
  gather_all_samples_fq
  merge_total_reads
  rm -rf "${TMP_DIR}/chunks" "$META_FOR_CHUNK"
  debug_kv "result.total_reads_rows" "$(wc -l <"$TOTAL_READS")"
}

split_fq_by_barcode() {
  chunk_all_input_fq
  gather_chunked_fq
}

demux_fastq_path() {
  local sample_name="$1"

  if [[ "$SIMULATION_MODE" = 1 ]]; then
    printf "%s/%s_all.fastq.gz\n" "$DEMUX_FASTQ_DIR" "$sample_name"
    return 0
  fi

  printf "%s/%s.fastq.gz\n" "$DEMUX_FASTQ_DIR" "$sample_name"
}

mapping_output_stem() {
  local sample_name="$1"
  local read_limit="$2"

  if [[ "$SIMULATION_MODE" = 1 ]]; then
    printf "%s_%s\n" "$sample_name" "$read_limit"
    return 0
  fi

  printf "%s\n" "$sample_name"
}

prepare_job_queue() {
  debug_log "STEP=prepare_job_queue"
  debug_kv "input.formatted_sample_sheet" "$FORMATTED_SAMPLE_SHEET"
  debug_kv "input.total_reads" "$TOTAL_READS"
  debug_kv "input.awk" "$PREPARE_JOB_QUEUE_AWK"
  debug_kv "output.job_queue" "$JOB_QUEUE"
  awk -F"\t" -v OFS="\t" -v sim_mode="$SIMULATION_MODE" \
    -v total_reads="$TOTAL_READS" -f "$PREPARE_JOB_QUEUE_AWK" \
    "$FORMATTED_SAMPLE_SHEET" >"$JOB_QUEUE"
  debug_kv "result.job_rows" "$(wc -l <"$JOB_QUEUE")"
}

read_mapping_job() {
  local job_id="$1"
  sed -n "${job_id}p" "$JOB_QUEUE" | cut -f2,4,6
}

build_mapping_paths() {
  local sample_name="$1"
  local read_limit="$2"
  local output_stem

  output_stem="$(mapping_output_stem "$sample_name" "$read_limit")"
  # scratch (fastq/bam/feature.tsv 等) は tmp/、診断ログは log/ に分離する。
  MAPPING_JOB_PREFIX="${TMP_DIR}/${output_stem}"
  MAPPING_JOB_LOG="${LOG_DIR}/${output_stem}.mapping_tools.log"
  MAPPING_JOB_SUMMARY="${LOG_DIR}/${output_stem}.mapping.log"
  MAPPING_JOB_FC_SUMMARY="${LOG_DIR}/${output_stem}.feature.tsv.summary"
  OUTPUT_BAM="${BAM_DIR}/${output_stem}.bam"
  INPUT_FASTQ="$(demux_fastq_path "$sample_name")"
}

reset_mapping_job_log() {
  : >"$MAPPING_JOB_LOG"
}

show_mapping_job_log_tail() {
  local sample_name="$1"
  local read_limit="$2"

  handle_error "Mapping failed for ${sample_name} (read count: ${read_limit})."
  handle_error "See log: ${MAPPING_JOB_LOG}"
  tail -n 20 "$MAPPING_JOB_LOG" >&2
}

prepare_mapping_input_fastq() {
  local read_limit="$1"

  if [[ "$read_limit" != "all" ]]; then
    seqkit head -n "$read_limit" "$INPUT_FASTQ" -o "${MAPPING_JOB_PREFIX}.fastq.gz"
    INPUT_FASTQ="${MAPPING_JOB_PREFIX}.fastq.gz"
  fi
}

run_trim_galore() {
  trim_galore -j "$THREADS" -o "$TMP_DIR" -a "GATCGTCGGACT" \
    --no_report_file --suppress_warn "$INPUT_FASTQ" >>"$MAPPING_JOB_LOG" 2>&1
}

cleanup_subsampled_fastq() {
  local read_limit="$1"

  if [[ "$read_limit" != "all" ]]; then
    rm -f "$INPUT_FASTQ"
  fi
}

# 実行時パラメータを summary/analysis_params.tsv に key<TAB>value で記録する。
# HTML レポートがこれを読み、Materials & Methods に実際のコマンドを表示する。
write_analysis_params() {
  local f="${SUMMARY_DIR}/analysis_params.tsv"
  {
    printf "key\tvalue\n"
    printf "hisat2_very_sensitive\t%s\n" "$HISAT2_VERY_SENSITIVE"
    printf "hisat2_score_min\t%s\n" "${HISAT2_SCORE_MIN:-}"
    printf "trim_adapter\t%s\n" "GATCGTCGGACT"
    printf "featurecounts_strand\t%s\n" "1"
  } >"$f"
}

run_hisat2_alignment() {
  local genome_name="$1"
  local hisat2_args=()
  local hisat2_index_prefix

  hisat2_index_prefix="$(pic_get_hisat2_index_prefix "$genome_name")"

  if [[ "$HISAT2_VERY_SENSITIVE" = 1 ]]; then
    hisat2_args+=(--very-sensitive)
  fi

  # --score-min のみを緩める運用 (--very-sensitive より高速で効果は同等)。
  # --very-sensitive と併用された場合は、後から渡すこちらが優先される。
  if [[ -n "${HISAT2_SCORE_MIN:-}" ]]; then
    hisat2_args+=(--score-min "$HISAT2_SCORE_MIN")
  fi

  {
    hisat2 -p "$THREADS" \
      -x "$hisat2_index_prefix" \
      "${hisat2_args[@]}" \
      -U "${MAPPING_JOB_PREFIX}_trimmed.fq.gz" \
      --summary-file "${MAPPING_JOB_SUMMARY}" |
      samtools view -@ "$THREADS" -huS |
      samtools sort -@ "$THREADS" -o "${MAPPING_JOB_PREFIX}.bam"
  } >>"$MAPPING_JOB_LOG" 2>&1
}

run_feature_counts() {
  local genome_name="$1"
  local gtf_path

  gtf_path="$(pic_get_gtf_path "$genome_name")"

  # featureCounts は summary を -o の隣 (tmp/) に書くため、生成後に log/ へ移す
  # (count 本体 .feature.tsv は scratch なので tmp/ のまま cleanup で破棄)。
  featureCounts -T "$THREADS" -s1 -a "$gtf_path" \
    -o "${MAPPING_JOB_PREFIX}.feature.tsv" -R BAM "${MAPPING_JOB_PREFIX}.bam" \
    >>"$MAPPING_JOB_LOG" 2>&1
  mv "${MAPPING_JOB_PREFIX}.feature.tsv.summary" "${MAPPING_JOB_FC_SUMMARY}"
}

finalize_mapping_bam() {
  samtools sort -@ "$THREADS" -o "$OUTPUT_BAM" \
    "${MAPPING_JOB_PREFIX}.bam.featureCounts.bam" >>"$MAPPING_JOB_LOG" 2>&1
  samtools index -@ "$THREADS" "$OUTPUT_BAM" >>"$MAPPING_JOB_LOG" 2>&1
}

cleanup_mapping_job_files() {
  rm -f "${MAPPING_JOB_PREFIX}"{_trimmed.fq.gz,.bam.featureCounts.bam,.bam,.feature.tsv}
}

run_single_mapping_job() {
  local genome_name="$1"
  local sample_name="$2"
  local read_limit="$3"

  build_mapping_paths "$sample_name" "$read_limit"
  reset_mapping_job_log
  log_info "Mapping ${sample_name} (read count: ${read_limit})"
  debug_log "STEP=run_single_mapping_job sample=${sample_name}"
  debug_kv "arg.genome" "$genome_name"
  debug_kv "arg.read_limit" "$read_limit"
  debug_kv "input.fastq" "$INPUT_FASTQ"
  debug_kv "output.bam" "$OUTPUT_BAM"
  debug_kv "output.job_log" "$MAPPING_JOB_LOG"

  prepare_mapping_input_fastq "$read_limit" || {
    show_mapping_job_log_tail "$sample_name" "$read_limit"
    return 1
  }

  run_trim_galore || {
    cleanup_subsampled_fastq "$read_limit"
    show_mapping_job_log_tail "$sample_name" "$read_limit"
    return 1
  }

  cleanup_subsampled_fastq "$read_limit"
  run_hisat2_alignment "$genome_name" || {
    show_mapping_job_log_tail "$sample_name" "$read_limit"
    return 1
  }

  run_feature_counts "$genome_name" || {
    show_mapping_job_log_tail "$sample_name" "$read_limit"
    return 1
  }

  finalize_mapping_bam || {
    show_mapping_job_log_tail "$sample_name" "$read_limit"
    return 1
  }

  cleanup_mapping_job_files
  debug_log "RESULT=run_single_mapping_job sample=${sample_name} status=ok"
}

exec_mapping_jobs() {
  local job_count job_id genome_name sample_name read_limit

  job_count="$(wc -l <"$JOB_QUEUE")"
  debug_log "STEP=exec_mapping_jobs"
  debug_kv "input.job_queue" "$JOB_QUEUE"
  debug_kv "arg.job_count" "$job_count"
  for job_id in $(seq "$job_count"); do
    read -r genome_name sample_name read_limit < <(read_mapping_job "$job_id")
    run_single_mapping_job "$genome_name" "$sample_name" "$read_limit" || return 1
  done
  debug_log "RESULT=exec_mapping_jobs status=ok"
}

map_and_count_feature() {
  prepare_job_queue
  exec_mapping_jobs
}

run_umi_count_for_bam() {
  local input_bam="$1"
  local output_table="$2"

  umi_tools count --method=unique --log2stderr \
    --per-gene --per-cell --gene-tag=XT \
    -I "$input_bam" -S "$output_table" >/dev/null 2>&1
}

umi_count() {
  local bam_path output_table

  while IFS= read -r bam_path; do
    wait_for_available_slot "$THREADS" || return 1
    output_table="${COUNTS_DIR}/$(basename "${bam_path%.bam}.txt")"
    run_umi_count_for_bam "$bam_path" "$output_table" &
  done < <(find "$BAM_DIR" -name "*.bam")

  wait
}

generate_deftables() {
  debug_log "STEP=generate_deftables"
  debug_kv "input.sample_sheet" "$SAMPLE_SHEET"
  debug_kv "input.counts_dir" "$COUNTS_DIR"
  debug_kv "input.awk" "$GENERATE_DEFTABLES_AWK"
  debug_kv "output.deftable_dir" "$DEFTABLE_DIR"
  mkdir -p "$DEFTABLE_DIR"

  tr -d "\015" <"$SAMPLE_SHEET" |
    awk -F"\t" -v OFS="\t" \
      -v counts_dir="$COUNTS_DIR" \
      -v deftable_dir="$DEFTABLE_DIR" \
      -v run_name="$RUN_NAME" \
      -v sim_mode="$SIMULATION_MODE" \
      -f "$GENERATE_DEFTABLES_AWK"
  debug_log "RESULT=generate_deftables status=ok"
}

mapping_sum() {
  debug_log "STEP=mapping_sum"
  debug_kv "input.total_reads" "$TOTAL_READS"
  debug_kv "input.counts_dir" "$COUNTS_DIR"
  debug_kv "input.job_queue" "$JOB_QUEUE"
  debug_kv "input.awk" "$MAP_SUMMARY_AWK"
  debug_kv "output.map_summary" "$MAP_SUMMARY"
  {
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "sample" "group" "barcode" "total" "trimmed" "unmapped" "multimapping" "nofeatures" \
      "ambiguity" "assigned" "umis" "genes" "umis/genes" "assigned/umis"

    awk -F"\t" -v OFS="\t" \
      -v log_dir="$LOG_DIR" -v count_dir="$COUNTS_DIR" \
      -v job_queue="$JOB_QUEUE" -f "$MAP_SUMMARY_AWK" \
      "$TOTAL_READS" | sort -k1,1 -k4,4n
  } >"$MAP_SUMMARY"
  debug_kv "result.map_summary_rows" "$(wc -l <"$MAP_SUMMARY")"

  if [[ "$SIMULATION_MODE" = 1 ]]; then
    /usr/local/bin/Rscript "$PLOT_SIM_SCRIPT" \
      "$MAP_SUMMARY" "$OUTPUT_DIR" >/dev/null 2>&1
  fi
}

chrom_sizes_file() {
  local genome_name="$1"
  pic_get_chromsize_path "$genome_name"
}

has_chrom_sizes_for_genome() {
  local genome_name="$1"
  [[ -f "$(chrom_sizes_file "$genome_name")" ]]
}

list_missing_chrom_sizes_genomes() {
  local genome_name

  cut -f2 "$JOB_QUEUE" | sort -u | while IFS= read -r genome_name; do
    if ! has_chrom_sizes_for_genome "$genome_name"; then
      printf "%s\n" "$genome_name"
    fi
  done
}

can_generate_bigwig() {
  [[ -f "$JOB_QUEUE" ]] || return 1
  [[ -z "$(list_missing_chrom_sizes_genomes)" ]]
}

run_umi_dedup() {
  local input_bam="$1"
  local output_bam="$2"

  umi_tools dedup --method=unique --log2stderr \
    --per-gene --per-cell --gene-tag=XT \
    -I "$input_bam" -S "$output_bam" >/dev/null 2>&1
}

index_umi_bam() {
  local umi_bam="$1"
  samtools index "$umi_bam" >/dev/null 2>&1
}

convert_to_umi() {
  local source_bam input_bam output_bam source_pattern

  if [[ "$SIMULATION_MODE" = 1 ]]; then
    source_pattern="*_all.bam"
  else
    source_pattern="*.bam"
  fi

  while IFS= read -r source_bam; do
    wait_for_available_slot "$THREADS" || return 1
    input_bam="$source_bam"
    if [[ "$SIMULATION_MODE" = 1 ]]; then
      output_bam="${UMI_BAM_DIR}/$(basename "${source_bam/_all.bam/_umi.bam}")"
    else
      output_bam="${UMI_BAM_DIR}/$(basename "${source_bam%.bam}")_umi.bam"
    fi
    run_umi_dedup "$input_bam" "$output_bam" &
  done < <(find "$BAM_DIR" -name "$source_pattern")
  wait || return 1

  while IFS= read -r output_bam; do
    wait_for_available_slot "$THREADS" || return 1
    index_umi_bam "$output_bam" &
  done < <(find "$UMI_BAM_DIR" -name "*_umi.bam")

  wait
}

single_bam_to_bw() {
  local job_id="$1"
  local umi_bam count bn fn genome sample

  read -r genome sample < <(
    sed -n "${job_id}p" "$JOB_QUEUE" | cut -f2,4
  )

  umi_bam="${UMI_BAM_DIR}/${sample}_umi.bam"
  count=$(samtools view -c "$umi_bam")
  bn="$(basename "$umi_bam")"
  fn=${bn/umi.bam/umi.cpm.bw}

  bedtools genomecov -split -bg -ibam "$umi_bam" | awk -v count="$count" '{
    $4 = $4 * 1000000 / count
    print
  }' | sort -k1,1 -k2,2n | sed 's/chrMT/chrM/g' >"$TMP_DIR/$fn.bg"

  bedGraphToBigWig "$TMP_DIR/$fn.bg" "$(chrom_sizes_file "$genome")" "$BIGWIG_DIR/$fn"
  rm -f "$TMP_DIR/$fn.bg"
}

bam_to_bw() {
  local job_queue

  mkdir -p "$TMP_DIR" "$BIGWIG_DIR"
  job_queue="$(wc -l <"$JOB_QUEUE")"
  run_numbered_jobs_in_parallel "$THREADS" "$job_queue" single_bam_to_bw
}

bam_to_cpm_bw() {
  convert_to_umi
  bam_to_bw
}

run_demux_command() {
  log_info "Validating inputs"
  validate_inputs
  log_info "Preparing demux workspace"
  prepare_workspace_for_demux
  log_info "Formatting sample sheet"
  format_sample_sheet
  if [[ "$SKIP_DEMUX" = 0 ]]; then
    log_info "Running demux by barcode"
    split_fq_by_barcode
  fi
  rm -rf "$TMP_DIR"
  debug_log "RESULT=run_demux_command status=ok"
}

run_primary_command() {
  debug_log "STEP=run_primary_command"
  debug_kv "arg.output_dir" "$OUTPUT_DIR"
  debug_kv "arg.sample_sheet" "$SAMPLE_SHEET"
  debug_kv "arg.run_name" "$RUN_NAME"
  debug_kv "arg.threads" "$THREADS"
  debug_kv "arg.skip_demux" "$SKIP_DEMUX"
  debug_kv "arg.input_fastq_dir" "$DEMUX_FASTQ_DIR"
  log_info "Validating inputs"
  validate_inputs
  validate_sample_sheet_genome_registration
  log_info "Preparing workspace"
  prepare_workspace_for_run
  # 実行時パラメータを記録する。HTML レポートの Materials & Methods は
  # このファイルを読んで、実際に使ったコマンドを表示する。
  write_analysis_params
  # scratch の tmp/ は関数終了時 (成功・失敗どちらでも) に掃除する。診断ログは
  # log/ に直接書かれるため、失敗しても log/ は残る。
  trap 'safe_rm_rf "$TMP_DIR"' RETURN
  log_info "Formatting sample sheet"
  format_sample_sheet

  if [[ "$SKIP_DEMUX" = 0 ]]; then
    log_info "Running demux by barcode"
    split_fq_by_barcode
  fi

  log_info "Running mapping and counting"
  map_and_count_feature
  log_info "Running UMI counting"
  umi_count
  log_info "Generating deftables and mapping summary"
  generate_deftables
  mapping_sum

  if can_generate_bigwig; then
    log_info "Generating BigWig files"
    bam_to_cpm_bw
  else
    handle_error "Skipping BigWig generation because one or more chrom size files are missing."
    handle_error "Run mapping again after preparing the required chrom size files."
  fi

  log_info "Finalizing outputs"
  finalize_outputs
  debug_log "RESULT=run_primary_command status=ok"
}
