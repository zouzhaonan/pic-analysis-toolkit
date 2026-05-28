#!/usr/bin/env bash

# 役割:
#   pic が使う genome ごとの参照ファイル登録を管理する。
# 入力:
#   genome 名と map TSV。
# 出力:
#   参照先 path の取得、登録一覧表示、map TSV 更新。

pic_genome_map_header() {
  printf "genome\tfasta_source\tgtf_source\n"
}

pic_effective_genome_map_file() {
  if [[ -n "${PIC_GENOME_MAP_FILE:-}" ]]; then
    printf "%s\n" "$PIC_GENOME_MAP_FILE"
  else
    printf "%s\n" "$(pic_resolve_lib_dir)/register/genome_map.tsv"
  fi
}

ensure_pic_genome_map_file() {
  local map_file
  local header_line
  local awk_dir
  awk_dir="${PIC_ROOT}/scripts/awk"
  map_file="$(pic_effective_genome_map_file)"
  mkdir -p "$(dirname "$map_file")"
  if [[ ! -f "$map_file" ]]; then
    pic_genome_map_header > "$map_file"
    return 0
  fi

  header_line="$(head -n 1 "$map_file" | tr -d '\r')"
  if [[ "$header_line" == *"\\t"* ]]; then
    awk -f "${awk_dir}/genome_map_fix_escaped_tabs.awk" \
      "$map_file" > "${map_file}.tmp"
    mv "${map_file}.tmp" "$map_file"
    header_line="$(head -n 1 "$map_file" | tr -d '\r')"
  fi

  if [[ "$header_line" != $'genome\tfasta_source\tgtf_source' ]]; then
    awk -f "${awk_dir}/genome_map_migrate_columns.awk" \
      "$map_file" > "${map_file}.tmp"
    mv "${map_file}.tmp" "$map_file"
  fi
}

read_pic_map_value() {
  local genome_name="$1"
  local column_name="$2"
  local map_file
  map_file="$(pic_effective_genome_map_file)"

  [[ -f "$map_file" ]] || return 0

  awk -F '\t' -v target="$genome_name" -v colname="$column_name" '
    BEGIN { gsub(/\r/, "", target) }
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == colname) {
          col = i
          break
        }
      }
      next
    }
    {
      g = $1
      gsub(/\r/, "", g)
    }
    col > 0 && g == target {
      v = $col
      gsub(/\r/, "", v)
      print v
      exit
    }
  ' "$map_file"
}

register_pic_genome_entry() {
  local genome_name="$1"
  local fasta_source="$2"
  local gtf_source="$3"
  local map_file
  local awk_dir
  awk_dir="${PIC_ROOT}/scripts/awk"
  map_file="$(pic_effective_genome_map_file)"

  ensure_pic_genome_map_file

  awk \
    -v genome="$genome_name" \
    -v fasta_src="$fasta_source" \
    -v gtf_src="$gtf_source" \
    -f "${awk_dir}/genome_map_upsert_entry.awk" \
    "$map_file" > "${map_file}.tmp"

  mv "${map_file}.tmp" "$map_file"
}

pic_get_hisat2_index_prefix() {
  local genome_name="$1"
  printf "%s\n" "$(pic_resolve_lib_dir)/hisat2_index/${genome_name}/genome"
}

pic_get_gtf_path() {
  local genome_name="$1"
  local mapped
  mapped="$(pic_resolve_lib_dir)/gtf/${genome_name}.gtf"
  if [[ -f "${mapped}.gz" ]]; then
    printf "%s\n" "${mapped}.gz"
    return 0
  fi
  printf "%s\n" "$mapped"
}

pic_get_chromsize_path() {
  local genome_name="$1"
  printf "%s\n" "$(pic_resolve_lib_dir)/chrom_size/${genome_name}.chrom.sizes"
}

pic_is_genome_registered() {
  local genome_name="$1"
  local map_file
  map_file="$(pic_effective_genome_map_file)"
  [[ -f "$map_file" ]] || return 1
  awk -F '\t' -v target="$genome_name" '
    BEGIN { gsub(/\r/, "", target) }
    NR>1 {
      g = $1
      gsub(/\r/, "", g)
      if (g == target) { found=1; exit }
    }
    END{exit(found?0:1)}
  ' "$map_file"
}

pic_validate_registered_genome() {
  local genome_name="$1"
  local idx gtf chrom
  local map_file
  map_file="$(pic_effective_genome_map_file)"

  idx="$(pic_get_hisat2_index_prefix "$genome_name")"
  gtf="$(pic_get_gtf_path "$genome_name")"
  chrom="$(pic_get_chromsize_path "$genome_name")"

  if ! pic_is_genome_registered "$genome_name"; then
    handle_error "Genome '${genome_name}' is not registered in ${map_file}."
    handle_error "Register it with: pic build-genome --genome ${genome_name} --fasta <path> --gtf <path>"
    return 1
  fi

  if [[ ! -f "${idx}.1.ht2" && ! -f "${idx}.1.ht2l" ]]; then
    handle_error "HISAT2 index was not found for '${genome_name}': ${idx}.1.ht2"
    return 1
  fi

  if [[ ! -f "$gtf" ]]; then
    handle_error "GTF was not found for '${genome_name}': ${gtf}"
    return 1
  fi

  if [[ ! -f "$chrom" ]]; then
    handle_error "chrom size file was not found for '${genome_name}': ${chrom}"
    return 1
  fi

  return 0
}

print_registered_pic_genomes() {
  local genomes
  local map_file
  map_file="$(pic_effective_genome_map_file)"
  echo "Registered pic genomes:"

  if [[ ! -f "$map_file" ]]; then
    echo "  (none)"
    return 0
  fi

  genomes="$(awk -F '\t' 'NR > 1 { g=$1; gsub(/\r/, "", g); if (g != "") print g }' "$map_file" | sort -u)"
  if [[ -z "$genomes" ]]; then
    echo "  (none)"
    return 0
  fi

  while IFS= read -r genome_name; do
    [[ -z "$genome_name" ]] && continue
    echo "  - ${genome_name}"
  done <<< "$genomes"
}
