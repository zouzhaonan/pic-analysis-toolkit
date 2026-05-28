#!/usr/bin/env bash

# 役割:
#   build-genome サブコマンド本体。FASTA から HISAT2 index を作り、
#   genome 参照情報を register する。
# 入力:
#   genome 名、fasta/gtf の path と source 情報。
# 出力:
#   HISAT2 index と register/genome_map.tsv。

run_build_genome_command() {
  local genome_name="$PREPARE_GENOME"
  local fasta_path="$PREPARE_FASTA"
  local gtf_path="$PREPARE_GTF"
  local fasta_source="$PREPARE_FASTA_SOURCE"
  local gtf_source="$PREPARE_GTF_SOURCE"
  local prepare_threads="$PREPARE_THREADS"
  local gtf_copy_path chromsize_copy_path

  ensure_pic_genome_map_file

  local index_dir index_prefix
  index_dir="${LIB_DIR}/hisat2_index/${genome_name}"
  index_prefix="${index_dir}/genome"
  gtf_copy_path="${LIB_DIR}/gtf/${genome_name}.gtf"
  if [[ "$gtf_path" == *.gz ]]; then
    gtf_copy_path="${gtf_copy_path}.gz"
  fi
  chromsize_copy_path="${LIB_DIR}/chrom_size/${genome_name}.chrom.sizes"

  mkdir -p "$index_dir" "$(dirname "$gtf_copy_path")" "$(dirname "$chromsize_copy_path")"

  log_info "Building HISAT2 index for ${genome_name}"
  hisat2-build -q -p "$prepare_threads" "$fasta_path" "$index_prefix"

  log_info "Copying GTF to ${gtf_copy_path}"
  cp -f "$gtf_path" "$gtf_copy_path"
  log_info "Generating chrom.sizes with faSize: ${chromsize_copy_path}"
  faSize -detailed -tab "$fasta_path" > "$chromsize_copy_path"

  register_pic_genome_entry \
    "$genome_name" \
    "$fasta_source" \
    "$gtf_source"

  log_info "Registered genome '${genome_name}' in ${PIC_GENOME_MAP_FILE}"
}
