#!/usr/bin/env bash

chipatlas_url_fixed="https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/"
chipatlas_poll_interval_fixed="60"
chipatlas_sbatch_options_fixed="-p epyc -t 180"
chipatlas_supported_genomes="hg38 mm10 rn6 dm6 ce11 sacCer3"

count_group_replicates() {
  local deftable_tsv="$1"
  local group_name="$2"

  awk -F '\t' -v target_group="$group_name" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "sample") sample_col = i
        if ($i == "group") group_col = i
      }
      next
    }
    {
      gsub(/\r$/, "", $group_col)
      gsub(/\r$/, "", $sample_col)
      if ($group_col == target_group && $sample_col != "") {
        seen[$sample_col] = 1
      }
    }
    END {
      count = 0
      for (sample_name in seen) count++
      print count
    }
  ' "$deftable_tsv"
}

count_nonempty_lines() {
  local file_path="$1"
  awk 'NF > 0 {count++} END {print count + 0}' "$file_path"
}

extract_chipatlas_job_url() {
  local response="$1"
  local job_url=""

  # 1) full URL
  job_url="$(printf "%s" "$response" | grep -Eo 'https?://[^[:space:]]*wabi_chipatlas_[^[:space:]]*' | head -n1 || true)"
  if [[ -n "$job_url" ]]; then
    printf "%s\n" "$job_url"
    return 0
  fi

  # 2) path-only URL
  local job_path=""
  job_path="$(printf "%s" "$response" | grep -Eo '/wabi/chipatlas/wabi_chipatlas_[^[:space:]]*' | head -n1 || true)"
  if [[ -n "$job_path" ]]; then
    printf "https://dtn1.ddbj.nig.ac.jp%s\n" "$job_path"
    return 0
  fi

  # 3) bare job id
  local job_id=""
  job_id="$(printf "%s" "$response" | grep -Eo 'wabi_chipatlas_[0-9A-Za-z._-]+' | head -n1 || true)"
  if [[ -n "$job_id" ]]; then
    printf "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/%s\n" "$job_id"
    return 0
  fi

  return 1
}

run_chipatlas_subcommand() {
  local deftable_tsv=""
  local umi_count_csv=""
  local deg_list_dir=""
  local genome=""
  local out_dir=""
  local antigen_class="TFs and others"
  local cell_class="All cell types"
  local threshold="100"
  local distance_up="5000"
  local distance_down="5000"

  while [[ $# -gt 0 ]]; do
    case "$1" in
    --help)
      show_help_for_subcommand "chipatlas"
      return 0
      ;;
    --deftable)
      deftable_tsv="$2"
      shift 2
      ;;
    --umi)
      umi_count_csv="$2"
      shift 2
      ;;
    --deg-list-dir)
      deg_list_dir="$2"
      shift 2
      ;;
    --genome)
      genome="$2"
      shift 2
      ;;
    --out-dir)
      out_dir="$2"
      shift 2
      ;;
    --antigen-class)
      antigen_class="$2"
      shift 2
      ;;
    --cell-class)
      cell_class="$2"
      shift 2
      ;;
    --threshold)
      threshold="$2"
      shift 2
      ;;
    --distance-up)
      distance_up="$2"
      shift 2
      ;;
    --distance-down)
      distance_down="$2"
      shift 2
      ;;
    *)
      handle_error "Unknown option for chipatlas: $1"
      return 1
      ;;
    esac
  done

  if [[ -z "$deftable_tsv" || -z "$umi_count_csv" || -z "$deg_list_dir" || -z "$genome" || -z "$out_dir" ]]; then
    handle_error "--deftable, --umi, --deg-list-dir, --genome, --out-dir are required."
    return 1
  fi

  if ! printf "%s\n" $chipatlas_supported_genomes | grep -Fxq "$genome"; then
    handle_error "Unsupported genome for ChIP-Atlas: ${genome}"
    handle_error "Supported genomes: hg38, mm10, rn6, dm6, ce11, sacCer3"
    return 1
  fi

  deftable_tsv="$(resolve_existing_path "$deftable_tsv")"
  umi_count_csv="$(resolve_existing_path "$umi_count_csv")"
  deg_list_dir="$(resolve_existing_path "$deg_list_dir")"
  out_dir="$(resolve_existing_path "$out_dir")"

  if [[ ! -f "$deftable_tsv" ]]; then
    handle_error "deftable TSV not found: $deftable_tsv"
    return 1
  fi
  if [[ ! -f "$umi_count_csv" ]]; then
    handle_error "UMI count CSV not found: $umi_count_csv"
    return 1
  fi
  if [[ ! -d "$deg_list_dir" ]]; then
    handle_error "DEGList directory not found: $deg_list_dir"
    return 1
  fi

  mkdir -p "$out_dir"

  debug_log "STEP=chipatlas"
  debug_kv "input.deftable" "$deftable_tsv"
  debug_kv "input.umi" "$umi_count_csv"
  debug_kv "input.deg_list_dir" "$deg_list_dir"
  debug_kv "arg.genome" "$genome"
  debug_kv "output.out_dir" "$out_dir"

  local pair_index=0
  local input_root="$out_dir/input"
  local result_root="$out_dir/result"
  local result_ea_root="$result_root/ea"
  local result_page_root="$result_root/page"
  local run_root="$out_dir/run"
  local jobs_tsv="$run_root/chipatlas_jobs.tsv"
  mkdir -p "$input_root" "$result_ea_root" "$result_page_root" "$run_root"
  : >"$jobs_tsv"
  local -a prefixes=()
  local file prefix

  while IFS= read -r file; do
    prefix="$(basename "$file")"
    prefix="${prefix%%__*}"
    prefixes+=("$prefix")
  done < <(find "$deg_list_dir" -maxdepth 1 -type f -name '*__*.txt' | sort)

  if [[ ${#prefixes[@]} -eq 0 ]]; then
    handle_error "No DEG list files were found in: $deg_list_dir"
    return 1
  fi

  local uniq_prefixes
  uniq_prefixes="$(printf "%s\n" "${prefixes[@]}" | sort -u)"

  while IFS= read -r prefix; do
    [[ -z "$prefix" ]] && continue
    pair_index=$((pair_index + 1))

    local -a pair_files=()
    while IFS= read -r file; do
      pair_files+=("$file")
    done < <(find "$deg_list_dir" -maxdepth 1 -type f -name "${prefix}__*.txt" | sort)

    if [[ ${#pair_files[@]} -ne 2 ]]; then
      handle_error "Expected 2 DEG lists for contrast ${prefix}, found ${#pair_files[@]}"
      continue
    fi

    local file_a="${pair_files[0]}"
    local file_b="${pair_files[1]}"
    local group_a group_b
    local group_a_reps group_b_reps

    group_a="$(basename "$file_a")"
    group_a="${group_a%.txt}"
    group_a="${group_a##*__}"
    group_b="$(basename "$file_b")"
    group_b="${group_b%.txt}"
    group_b="${group_b##*__}"
    group_a_reps="$(count_group_replicates "$deftable_tsv" "$group_a")"
    group_b_reps="$(count_group_replicates "$deftable_tsv" "$group_b")"

    log_info "[ChIP-Atlas] (${pair_index}) Running ${prefix}: ${group_a} vs ${group_b}"
    debug_kv "chipatlas.${prefix}.replicates.${group_a}" "$group_a_reps"
    debug_kv "chipatlas.${prefix}.replicates.${group_b}" "$group_b_reps"

    local deg_count_a deg_count_b
    deg_count_a="$(count_nonempty_lines "$file_a")"
    deg_count_b="$(count_nonempty_lines "$file_b")"
    debug_kv "chipatlas.${prefix}.deg_count.${group_a}" "$deg_count_a"
    debug_kv "chipatlas.${prefix}.deg_count.${group_b}" "$deg_count_b"

    local submit_ea=1
    local submit_page=1
    local input_contrast_dir="$input_root/$prefix"
    local result_ea_dir="$result_ea_root/$prefix"
    local result_page_dir="$result_page_root/$prefix"

    if [[ "$deg_count_a" -eq 0 || "$deg_count_b" -eq 0 ]]; then
      submit_ea=0
      log_info "[ChIP-Atlas] Skipping EA for ${prefix}: empty DEG list detected (${group_a}=${deg_count_a}, ${group_b}=${deg_count_b})"
    fi

    if [[ "$group_a_reps" -lt 2 || "$group_b_reps" -lt 2 ]]; then
      submit_page=0
      log_info "[ChIP-Atlas] Skipping PAGE for ${prefix}: replicate count must be >= 2 in both groups (${group_a}=${group_a_reps}, ${group_b}=${group_b_reps})"
    fi

    if [[ "$submit_ea" -eq 0 && "$submit_page" -eq 0 ]]; then
      continue
    fi

    if [[ "$submit_ea" -eq 1 ]]; then
      mkdir -p "$input_contrast_dir" "$result_ea_dir"
      cp "$file_a" "$input_contrast_dir/${prefix}__${group_a}.txt"
      cp "$file_b" "$input_contrast_dir/${prefix}__${group_b}.txt"

      local ea_url
      ea_url="$(submit_chipatlas_ea \
        "$genome" "$antigen_class" "$cell_class" "$threshold" \
        "$distance_up" "$distance_down" \
        "$prefix" "$group_a" "$group_b" \
        "$input_contrast_dir/${prefix}__${group_a}.txt" "$input_contrast_dir/${prefix}__${group_b}.txt" \
        "$result_ea_dir")" || {
        handle_error "EA submit failed: ${prefix}"
        continue
      }
      printf "%s\t%s\t%s\t%s\n" "$prefix" "ea" "$ea_url" "${result_ea_dir}/${prefix}__chipatlas_ea" >>"$jobs_tsv"
    fi

    if [[ "$submit_page" -eq 1 ]]; then
      mkdir -p "$input_contrast_dir" "$result_page_dir"
      local page_input_csv="$input_contrast_dir/${prefix}__page_input.csv"
      build_chipatlas_page_input "$umi_count_csv" "$deftable_tsv" "$group_a" "$group_b" "$page_input_csv"

      local page_url
      page_url="$(submit_chipatlas_page \
        "$genome" "$antigen_class" "$cell_class" "$threshold" \
        "$distance_up" "$distance_down" \
        "$prefix" "$page_input_csv" "$result_page_dir")" || {
        handle_error "PAGE submit failed: ${prefix}"
        continue
      }
      printf "%s\t%s\t%s\t%s\n" "$prefix" "page" "$page_url" "${result_page_dir}/${prefix}__chipatlas_page" >>"$jobs_tsv"
    fi
  done <<<"$uniq_prefixes"

  run_chipatlas_jobs_in_parallel "$jobs_tsv" || return 1
  rm -rf "$run_root"
  log_info "[ChIP-Atlas] Completed. Output: $out_dir"
}

build_chipatlas_page_input() {
  local umi_count_csv="$1"
  local deftable_tsv="$2"
  local group_a="$3"
  local group_b="$4"
  local out_csv="$5"

  Rscript -e '
args <- commandArgs(trailingOnly=TRUE)
count_file <- args[[1]]
def_file <- args[[2]]
g1 <- args[[3]]
g2 <- args[[4]]
out_file <- args[[5]]

count_tbl <- read.csv(count_file, check.names = FALSE, stringsAsFactors = FALSE)
def_tbl <- read.delim(def_file, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("sample", "group") %in% colnames(def_tbl))) stop("deftable must include sample/group")
if (!("ens_gene" %in% colnames(count_tbl))) stop("UMI count CSV must include ens_gene")

samples <- unique(def_tbl$sample[def_tbl$group %in% c(g1, g2)])
if (length(samples) == 0) stop(sprintf("No samples for groups: %s, %s", g1, g2))

missing_cols <- setdiff(samples, colnames(count_tbl))
if (length(missing_cols) > 0) stop(sprintf("Missing sample columns in UMI count CSV: %s", paste(missing_cols, collapse=",")))

gene <- if ("ext_gene" %in% colnames(count_tbl)) count_tbl$ext_gene else count_tbl$ens_gene
gene[is.na(gene) | gene == "" | gene == "NA"] <- count_tbl$ens_gene[is.na(gene) | gene == "" | gene == "NA"]

out <- data.frame(gene = gene, count_tbl[, samples, drop = FALSE], check.names = FALSE)
out <- out[!is.na(out$gene) & out$gene != "", , drop = FALSE]
write.csv(out, out_file, row.names = FALSE, quote = FALSE)
' "$umi_count_csv" "$deftable_tsv" "$group_a" "$group_b" "$out_csv"
}

submit_chipatlas_ea() {
  local genome="$1"
  shift
  local antigen_class="$1"
  shift
  local cell_class="$1"
  shift
  local threshold="$1"
  shift
  local distance_up="$1"
  shift
  local distance_down="$1"
  shift
  local title="$1"
  shift
  local group_a="$1"
  shift
  local group_b="$1"
  shift
  local file_a="$1"
  shift
  local file_b="$1"
  shift
  local out_dir="$1"

  local response job_url
  local -a curl_args=(
    -sS -X POST
    -d "format=text"
    -d "result=www"
    -d "genome=${genome}"
    -d "antigenClass=${antigen_class}"
    -d "cellClass=${cell_class}"
    -d "threshold=${threshold}"
    -d "typeA=gene"
    --data-urlencode "bedAFile@${file_a}"
    -d "typeB=gene"
    --data-urlencode "bedBFile@${file_b}"
    -d "permTime=1"
    -d "title=${title}"
    -d "descriptionA=${group_a}"
    -d "descriptionB=${group_b}"
    -d "distanceUp=${distance_up}"
    -d "distanceDown=${distance_down}"
    -d "sbatchOptions=${chipatlas_sbatch_options_fixed}"
    "${chipatlas_url_fixed}"
  )
  debug_command curl "${curl_args[@]}"
  response="$(curl "${curl_args[@]}")"

  job_url="$(extract_chipatlas_job_url "$response" || true)"
  debug_kv "chipatlas.ea.response" "$(printf "%s" "$response" | tr '\n' ' ' | cut -c1-300)"
  debug_kv "chipatlas.ea.job_url" "$job_url"
  if [[ -z "$job_url" ]]; then
    handle_error "Failed to submit ChIP-Atlas EA job for ${title}"
    return 1
  fi
  printf "%s\n" "$job_url"
}

submit_chipatlas_page() {
  local genome="$1"
  shift
  local antigen_class="$1"
  shift
  local cell_class="$1"
  shift
  local threshold="$1"
  shift
  local distance_up="$1"
  shift
  local distance_down="$1"
  shift
  local title="$1"
  shift
  local page_input_csv="$1"
  shift
  local out_dir="$1"

  local response job_url
  local -a curl_args=(
    -sS -X POST
    -d "format=text"
    -d "result=www"
    -d "genome=${genome}"
    -d "antigenClass=${antigen_class}"
    -d "cellClass=${cell_class}"
    -d "threshold=${threshold}"
    -d "typeA=count"
    --data-urlencode "bedAFile@${page_input_csv}"
    -d "typeB=empty"
    -d "bedBFile=empty"
    -d "permTime=1"
    -d "title=${title}"
    -d "descriptionA=empty"
    -d "descriptionB=empty"
    -d "distanceUp=${distance_up}"
    -d "distanceDown=${distance_down}"
    -d "sbatchOptions=${chipatlas_sbatch_options_fixed}"
    "${chipatlas_url_fixed}"
  )
  debug_command curl "${curl_args[@]}"
  response="$(curl "${curl_args[@]}")"

  job_url="$(extract_chipatlas_job_url "$response" || true)"
  debug_kv "chipatlas.page.response" "$(printf "%s" "$response" | tr '\n' ' ' | cut -c1-300)"
  debug_kv "chipatlas.page.job_url" "$job_url"
  if [[ -z "$job_url" ]]; then
    handle_error "Failed to submit ChIP-Atlas PAGE job for ${title}"
    return 1
  fi
  printf "%s\n" "$job_url"
}

wait_and_download_chipatlas() {
  local job_url="$1"
  local out_prefix="$2"

  while :; do
    local status
    status="$(curl -s "$job_url" | awk '$1 == "status:" {print $2}')"

    if [[ "$status" == "finished" ]]; then
      curl -s "${job_url}?info=result&format=tsv" >"${out_prefix}.tsv"
      curl -s "${job_url}?info=result&format=html" >"${out_prefix}.html"
      curl -s "${job_url}?info=result&format=log" >"${out_prefix}.log"
      break
    fi

    if [[ "$status" == "error" || "$status" == "failed" ]]; then
      handle_error "ChIP-Atlas job failed: ${job_url}"
      return 1
    fi

    sleep "$chipatlas_poll_interval_fixed"
  done
}

run_chipatlas_jobs_in_parallel() {
  local jobs_tsv="$1"
  local line prefix job_type job_url out_prefix
  local -a pids=()
  local fail=0

  while IFS=$'\t' read -r prefix job_type job_url out_prefix; do
    [[ -z "${job_url:-}" ]] && continue
    debug_log "chipatlas.wait.start contrast=${prefix} type=${job_type} url=${job_url}"
    (
      wait_and_download_chipatlas "$job_url" "$out_prefix"
    ) &
    pids+=("$!")
  done <"$jobs_tsv"

  local pid
  for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
      fail=1
    fi
  done

  if [[ "$fail" -ne 0 ]]; then
    handle_error "One or more ChIP-Atlas jobs failed."
    return 1
  fi
}
