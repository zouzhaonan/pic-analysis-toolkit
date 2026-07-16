#!/usr/bin/env bash

# 役割:
#   xengsort を使った xenograft リード分類機能。
#     - pic build-xenograft-index : host/graft ゲノムから xengsort index を作成・登録
#     - pic xenograft             : demux 済み FASTQ を host/graft に分類し、
#                                   マッピングへ渡せる形に配置する
# 入力:
#   index 作成: host/graft の FASTA (DNA(+cDNA))、host/graft の pic genome 名。
#   分類:       登録済み index 名、demux 済み FASTQ ディレクトリ。
# 出力:
#   index 作成: <PIC_LIB>/xengsort_index/<name>/<name>.hash|.info と
#               register/xengsort_index_map.tsv への登録。
#   分類:       <out>/xengsort/<sample>-{host,graft,both,neither,ambiguous}.fq.gz、
#               <out>/classified/<fraction>/<sample>.fastq.gz、
#               xenograft_classify_summary__<run>.tsv (分類 QC)。
# 注記:
#   xengsort / ntcard は pic の conda パッケージに同梱され PATH 上にある前提。

# ---------------------------------------------------------------------------
# 依存コマンドの確認
# ---------------------------------------------------------------------------

xenograft_require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    handle_error "Required command not found: ${cmd}"
    handle_error "xengsort/ntcard は pic の conda パッケージに含まれます。環境を確認してください。"
    exit 1
  fi
}

# ---------------------------------------------------------------------------
# xengsort index レジストリ
#   register/xengsort_index_map.tsv
#   列: name  host_genome  graft_genome  k  n  host_source  graft_source
# ---------------------------------------------------------------------------

pic_xengsort_register_file() {
  printf "%s\n" "$(pic_resolve_lib_dir)/register/xengsort_index_map.tsv"
}

pic_xengsort_index_dir() {
  printf "%s\n" "$(pic_resolve_lib_dir)/xengsort_index/$1"
}

# xengsort classify --index に渡す prefix (拡張子なし)
pic_xengsort_index_prefix() {
  printf "%s\n" "$(pic_xengsort_index_dir "$1")/$1"
}

pic_xengsort_register_header() {
  printf "name\thost_genome\tgraft_genome\tk\tn\thost_source\tgraft_source\n"
}

ensure_xengsort_register_file() {
  local f
  f="$(pic_xengsort_register_file)"
  mkdir -p "$(dirname "$f")"
  if [[ ! -f "$f" ]]; then
    pic_xengsort_register_header > "$f"
  fi
}

# name をキーに upsert (既存行を除去して追記)
register_xengsort_index_entry() {
  local name="$1" host_genome="$2" graft_genome="$3" k="$4" n="$5" host_src="$6" graft_src="$7"
  local f tmp
  ensure_xengsort_register_file
  f="$(pic_xengsort_register_file)"
  tmp="${f}.tmp"

  awk -F '\t' -v name="$name" 'NR==1 || $1 != name' "$f" > "$tmp"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$name" "$host_genome" "$graft_genome" "$k" "$n" "$host_src" "$graft_src" >> "$tmp"
  mv "$tmp" "$f"
}

pic_xengsort_read_value() {
  local name="$1" col="$2" f
  f="$(pic_xengsort_register_file)"
  [[ -f "$f" ]] || return 0
  awk -F '\t' -v target="$name" -v colname="$col" '
    NR==1 { for (i=1;i<=NF;i++) if ($i==colname) c=i; next }
    $1==target && c>0 { print $c; exit }
  ' "$f"
}

pic_xengsort_is_registered() {
  local name="$1" f
  f="$(pic_xengsort_register_file)"
  [[ -f "$f" ]] || return 1
  awk -F '\t' -v t="$name" 'NR>1 && $1==t {found=1; exit} END{exit(found?0:1)}' "$f"
}

print_registered_xengsort_indexes() {
  local f
  f="$(pic_xengsort_register_file)"
  echo "Registered xengsort indexes:"
  if [[ ! -f "$f" ]]; then
    echo "  (none)"
    return 0
  fi
  local any=0
  while IFS=$'\t' read -r name host graft k n _ _; do
    [[ -z "$name" || "$name" == "name" ]] && continue
    echo "  - ${name} (host=${host}, graft=${graft}, k=${k}, n=${n})"
    any=1
  done < "$f"
  [[ "$any" = 1 ]] || echo "  (none)"
}

# ---------------------------------------------------------------------------
# pic build-xenograft-index
# ---------------------------------------------------------------------------

# FASTA 群 (.gz でも可) を 1 本の .fa.gz に結合する。
xenograft_combine_fastas() {
  local out_gz="$1"
  shift
  : > "$out_gz"
  local fa
  for fa in "$@"; do
    [[ -z "$fa" ]] && continue
    if [[ "$fa" == *.gz ]]; then
      cat "$fa" >> "$out_gz"
    else
      pigz -c "$fa" >> "$out_gz"
    fi
  done
}

# cDNA FASTA が無い生物種向けに、GTF + ゲノム FASTA から gffread で cDNA を生成する。
xenograft_cdna_from_gtf() {
  local genome_fa="$1" gtf="$2" out_cdna="$3" threads="$4" tmp_dir="$5" tag="$6"
  local g_fa="$genome_fa" g_gtf="$gtf"
  if [[ "$genome_fa" == *.gz ]]; then
    g_fa="${tmp_dir}/${tag}_genome.fa"
    pigz -dc -p "$threads" "$genome_fa" > "$g_fa"
  fi
  if [[ "$gtf" == *.gz ]]; then
    g_gtf="${tmp_dir}/${tag}.gtf"
    pigz -dc -p "$threads" "$gtf" > "$g_gtf"
  fi
  gffread "$g_gtf" -g "$g_fa" -w "$out_cdna" >&2
}

# ntcard で distinct k-mer 数 (F0) を推定する。
xenograft_estimate_n() {
  local k="$1" threads="$2" tmp_dir="$3"
  shift 3
  local plain=()
  local fa i=0
  for fa in "$@"; do
    i=$((i + 1))
    local p="${tmp_dir}/nt_${i}.fa"
    if [[ "$fa" == *.gz ]]; then
      pigz -dc -p "$threads" "$fa" > "$p"
    else
      cp -f "$fa" "$p"
    fi
    plain+=("$p")
  done
  ntcard -k "$k" -t "$threads" -p "${tmp_dir}/nt" "${plain[@]}" >&2
  local hist="${tmp_dir}/nt_k${k}.hist"
  awk '$1=="F0"{print $2; exit}' "$hist"
}

run_build_xenograft_index_subcommand() {
  if [[ $# -eq 0 ]]; then
    show_help_for_subcommand build-xenograft-index
    return 0
  fi

  local name="" host_genome="" graft_genome=""
  local host_fasta="" graft_fasta="" host_cdna="" graft_cdna="" host_gtf="" graft_gtf=""
  local threads="$PIC_DEFAULT_THREADS"
  local host_source="unspecified" graft_source="unspecified"
  local k=25 n=""   # k は xengsort 標準の 25 固定、n は ntcard で推定

  while [[ $# -gt 0 ]]; do
    case "$1" in
    --name) name="${2:-}"; shift 2 ;;
    --host-genome) host_genome="${2:-}"; shift 2 ;;
    --graft-genome) graft_genome="${2:-}"; shift 2 ;;
    --host-fasta) host_fasta="${2:-}"; shift 2 ;;
    --graft-fasta) graft_fasta="${2:-}"; shift 2 ;;
    --host-cdna) host_cdna="${2:-}"; shift 2 ;;
    --graft-cdna) graft_cdna="${2:-}"; shift 2 ;;
    --host-gtf) host_gtf="${2:-}"; shift 2 ;;
    --graft-gtf) graft_gtf="${2:-}"; shift 2 ;;
    --threads) threads="${2:-}"; shift 2 ;;
    --host-source) host_source="${2:-}"; shift 2 ;;
    --graft-source) graft_source="${2:-}"; shift 2 ;;
    --help) show_help_for_subcommand build-xenograft-index; return 0 ;;
    *) handle_error "Unknown option: $1"; exit 1 ;;
    esac
  done

  if [[ -z "$name" || -z "$host_genome" || -z "$graft_genome" || -z "$host_fasta" || -z "$graft_fasta" ]]; then
    handle_error "build-xenograft-index requires --name, --host-genome, --graft-genome, --host-fasta, --graft-fasta."
    exit 1
  fi
  local fa
  for fa in "$host_fasta" "$graft_fasta" \
    ${host_cdna:+"$host_cdna"} ${graft_cdna:+"$graft_cdna"} \
    ${host_gtf:+"$host_gtf"} ${graft_gtf:+"$graft_gtf"}; do
    if [[ ! -f "$fa" ]]; then
      handle_error "File not found: ${fa}"
      exit 1
    fi
  done

  xenograft_require_command xengsort
  xenograft_require_command pigz
  xenograft_require_command ntcard
  if [[ ( -z "$host_cdna" && -n "$host_gtf" ) || ( -z "$graft_cdna" && -n "$graft_gtf" ) ]]; then
    xenograft_require_command gffread
  fi

  local index_dir prefix tmp_dir
  index_dir="$(pic_xengsort_index_dir "$name")"
  prefix="$(pic_xengsort_index_prefix "$name")"
  mkdir -p "$index_dir"
  tmp_dir="$(mktemp -d)"
  # shellcheck disable=SC2064
  trap "rm -rf '$tmp_dir'" RETURN

  # cDNA が未指定で GTF がある場合は gffread で自動生成する。
  if [[ -z "$host_cdna" && -n "$host_gtf" ]]; then
    log_info "xenograft: host cDNA を GTF から生成します (gffread)"
    host_cdna="${tmp_dir}/host.cdna.fa"
    xenograft_cdna_from_gtf "$host_fasta" "$host_gtf" "$host_cdna" "$threads" "$tmp_dir" host
  fi
  if [[ -z "$graft_cdna" && -n "$graft_gtf" ]]; then
    log_info "xenograft: graft cDNA を GTF から生成します (gffread)"
    graft_cdna="${tmp_dir}/graft.cdna.fa"
    xenograft_cdna_from_gtf "$graft_fasta" "$graft_gtf" "$graft_cdna" "$threads" "$tmp_dir" graft
  fi

  log_info "xenograft: host/graft FASTA を結合します"
  local host_combined="${tmp_dir}/host.fa.gz" graft_combined="${tmp_dir}/graft.fa.gz"
  xenograft_combine_fastas "$host_combined" "$host_fasta" ${host_cdna:+"$host_cdna"}
  xenograft_combine_fastas "$graft_combined" "$graft_fasta" ${graft_cdna:+"$graft_cdna"}

  log_info "xenograft: ntcard で k-mer 数 (F0) を推定します (k=${k})"
  n="$(xenograft_estimate_n "$k" "$threads" "$tmp_dir" "$host_combined" "$graft_combined")"
  if [[ -z "$n" ]]; then
    handle_error "ntcard による F0 推定に失敗しました。入力 FASTA を確認してください。"
    exit 1
  fi
  log_info "xenograft: 推定 F0 (= -n) = ${n}"

  log_info "xenograft: xengsort index を作成します (name=${name}, host=-H, graft=-G)"
  xengsort index \
    --index "$prefix" \
    -H "$host_combined" \
    -G "$graft_combined" \
    -k "$k" \
    -n "$n" \
    --threads-read "$threads" \
    --threads-split "$threads" \
    -W "$threads"

  if [[ ! -f "${prefix}.hash" || ! -f "${prefix}.info" ]]; then
    handle_error "xengsort index の出力が見つかりません: ${prefix}.hash / .info"
    exit 1
  fi

  register_xengsort_index_entry "$name" "$host_genome" "$graft_genome" "$k" "$n" "$host_source" "$graft_source"
  log_info "xenograft: index '${name}' を登録しました ($(pic_xengsort_register_file))"
  log_info "  host_genome=${host_genome}  graft_genome=${graft_genome}"
}

# ---------------------------------------------------------------------------
# pic xenograft (classify)
# ---------------------------------------------------------------------------

# gzip FASTQ のリード数 (= 行数 / 4)
xenograft_count_reads() {
  local f="$1"
  if [[ ! -s "$f" ]]; then
    printf "0\n"
    return 0
  fi
  printf "%s\n" "$(( $(pigz -dc "$f" | wc -l) / 4 ))"
}

# 元の sample sheet から画分用 sample sheet を生成する (genome を差し替え、
# fastq_prefix/sample は分類 FASTQ 名 = 元の sample 列に合わせる)。
xenograft_write_fraction_sheet() {
  local orig="$1" genome="$2" out_sheet="$3"
  awk -F '\t' -v G="$genome" '
    BEGIN { OFS="\t"; print "fastq_prefix","genome","barcode","sample","group" }
    NR==1 && tolower($1)=="fastq_prefix" { next }
    NF>=5 { print $4, G, $3, $4, $5 }
  ' "$orig" > "$out_sheet"
}

run_xenograft_subcommand() {
  if [[ $# -eq 0 ]]; then
    show_help_for_subcommand xenograft
    return 0
  fi

  local index_name="" demux_dir="" out_dir="" fraction="both"
  local raw_fastq_dir="" sample_sheet=""
  local threads="$PIC_DEFAULT_THREADS" run_name=""

  while [[ $# -gt 0 ]]; do
    case "$1" in
    --index) index_name="${2:-}"; shift 2 ;;
    --demux-fastq-dir) demux_dir="${2:-}"; shift 2 ;;
    --raw-fastq-dir) raw_fastq_dir="${2:-}"; shift 2 ;;
    --sample-sheet) sample_sheet="${2:-}"; shift 2 ;;
    --out-dir) out_dir="${2:-}"; shift 2 ;;
    --fraction) fraction="${2:-}"; shift 2 ;;
    --threads) threads="${2:-}"; shift 2 ;;
    --run-name) run_name="${2:-}"; shift 2 ;;
    --help) show_help_for_subcommand xenograft; return 0 ;;
    *) handle_error "Unknown option: $1"; exit 1 ;;
    esac
  done

  if [[ -z "$index_name" || -z "$out_dir" ]]; then
    handle_error "xenograft requires --index and --out-dir."
    exit 1
  fi
  # 入力は (A) 既に demux 済みの --demux-fastq-dir、
  #        (B) --raw-fastq-dir + --sample-sheet (内部で demux) のどちらか。
  if [[ -z "$demux_dir" && ( -z "$raw_fastq_dir" || -z "$sample_sheet" ) ]]; then
    handle_error "xenograft には --demux-fastq-dir、または --raw-fastq-dir + --sample-sheet が必要です。"
    exit 1
  fi
  case "$fraction" in
  graft|host|both) ;;
  *) handle_error "--fraction は graft / host / both のいずれかです: ${fraction}"; exit 1 ;;
  esac
  if [[ -n "$sample_sheet" && ! -f "$sample_sheet" ]]; then
    handle_error "sample-sheet not found: ${sample_sheet}"
    exit 1
  fi

  xenograft_require_command xengsort
  xenograft_require_command pigz

  if ! pic_xengsort_is_registered "$index_name"; then
    handle_error "xengsort index '${index_name}' は未登録です。"
    handle_error "作成: pic build-xenograft-index --name ${index_name} ..."
    exit 1
  fi
  local prefix host_genome graft_genome
  prefix="$(pic_xengsort_index_prefix "$index_name")"
  host_genome="$(pic_xengsort_read_value "$index_name" host_genome)"
  graft_genome="$(pic_xengsort_read_value "$index_name" graft_genome)"
  if [[ ! -f "${prefix}.hash" || ! -f "${prefix}.info" ]]; then
    handle_error "index ファイルが見つかりません: ${prefix}.hash / .info"
    exit 1
  fi

  [[ -n "$run_name" ]] || run_name="$(basename "$out_dir")"
  mkdir -p "$out_dir"

  # モード B: raw FASTQ を内部で demux して per-sample FASTQ を得る。
  if [[ -z "$demux_dir" ]]; then
    demux_dir="${out_dir}/demux"
    log_info "xenograft: raw FASTQ を demux します -> ${demux_dir}"
    if ! "${PIC_ROOT}/bin/pic" demux \
      --sample-sheet "$sample_sheet" \
      --raw-fastq-dir "$raw_fastq_dir" \
      --out-dir "$demux_dir"; then
      handle_error "xenograft: demux に失敗しました。"
      exit 1
    fi
  fi
  if [[ ! -d "$demux_dir" ]]; then
    handle_error "demux-fastq-dir not found: ${demux_dir}"
    exit 1
  fi

  local xg_dir="${out_dir}/xengsort"
  mkdir -p "$xg_dir"
  local summary="${out_dir}/xenograft_classify_summary__${run_name}.tsv"
  printf "sample\thost\tgraft\tboth\tneither\tambiguous\ttotal\n" > "$summary"

  # 分類対象の FASTQ を収集 (.fastq.gz / .fq.gz)
  local fastqs=()
  local f
  for f in "$demux_dir"/*.fastq.gz "$demux_dir"/*.fq.gz; do
    [[ -e "$f" ]] && fastqs+=("$f")
  done
  if [[ "${#fastqs[@]}" -eq 0 ]]; then
    handle_error "FASTQ が見つかりません: ${demux_dir}/*.fastq.gz|*.fq.gz"
    exit 1
  fi

  log_info "xenograft: ${#fastqs[@]} 個の FASTQ を分類します (index=${index_name})"

  local sample base
  for f in "${fastqs[@]}"; do
    base="$(basename "$f")"
    sample="${base%.fastq.gz}"
    sample="${sample%.fq.gz}"
    log_info "xenograft: classify ${sample}"
    xengsort classify \
      --index "$prefix" \
      --fastq "$f" \
      -o "${xg_dir}/${sample}" \
      --mode count \
      -T "$threads"

    local c_host c_graft c_both c_neither c_ambi c_total
    c_host="$(xenograft_count_reads "${xg_dir}/${sample}-host.fq.gz")"
    c_graft="$(xenograft_count_reads "${xg_dir}/${sample}-graft.fq.gz")"
    c_both="$(xenograft_count_reads "${xg_dir}/${sample}-both.fq.gz")"
    c_neither="$(xenograft_count_reads "${xg_dir}/${sample}-neither.fq.gz")"
    c_ambi="$(xenograft_count_reads "${xg_dir}/${sample}-ambiguous.fq.gz")"
    c_total=$((c_host + c_graft + c_both + c_neither + c_ambi))
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$sample" "$c_host" "$c_graft" "$c_both" "$c_neither" "$c_ambi" "$c_total" >> "$summary"
  done

  # 画分ごとに、マッピングへ渡す FASTQ を配置する。
  # 規約: graft -> graft genome, host -> host genome。
  local fractions=()
  case "$fraction" in
  both) fractions=(graft host) ;;
  graft) fractions=(graft) ;;
  host) fractions=(host) ;;
  esac
  local fr genome_for
  for fr in "${fractions[@]}"; do
    local fr_dir="${out_dir}/classified/${fr}"
    mkdir -p "$fr_dir"
    for f in "${fastqs[@]}"; do
      base="$(basename "$f")"
      sample="${base%.fastq.gz}"
      sample="${sample%.fq.gz}"
      local src="${xg_dir}/${sample}-${fr}.fq.gz"
      [[ -f "$src" ]] || continue
      cp -f "$src" "${fr_dir}/${sample}.fastq.gz"
    done
    if [[ "$fr" == "graft" ]]; then genome_for="$graft_genome"; else genome_for="$host_genome"; fi
    log_info "xenograft: ${fr} 画分を配置しました -> ${fr_dir}  (genome=${genome_for})"
  done

  # sample sheet が渡されていれば、画分ごとの sample sheet を自動生成する。
  local graft_sheet="" host_sheet=""
  if [[ -n "$sample_sheet" ]]; then
    if [[ "$fraction" == "both" || "$fraction" == "graft" ]]; then
      graft_sheet="${out_dir}/sample_sheet_graft.tsv"
      xenograft_write_fraction_sheet "$sample_sheet" "$graft_genome" "$graft_sheet"
      log_info "xenograft: graft sample sheet を生成 -> ${graft_sheet} (genome=${graft_genome})"
    fi
    if [[ "$fraction" == "both" || "$fraction" == "host" ]]; then
      host_sheet="${out_dir}/sample_sheet_host.tsv"
      xenograft_write_fraction_sheet "$sample_sheet" "$host_genome" "$host_sheet"
      log_info "xenograft: host sample sheet を生成 -> ${host_sheet} (genome=${host_genome})"
    fi
  fi

  log_info "xenograft: 分類サマリ: ${summary}"
  log_info "xenograft: 次のステップ (画分ごとに pic all):"
  if [[ "$fraction" == "both" || "$fraction" == "graft" ]]; then
    log_info "  graft (genome=${graft_genome}): pic all --no-report --demux-fastq-dir ${out_dir}/classified/graft --sample-sheet ${graft_sheet:-<graft sheet>} --run-name ${run_name}_graft --out-dir ${out_dir}/graft"
  fi
  if [[ "$fraction" == "both" || "$fraction" == "host" ]]; then
    log_info "  host  (genome=${host_genome}): pic all --no-report --demux-fastq-dir ${out_dir}/classified/host  --sample-sheet ${host_sheet:-<host sheet>}  --run-name ${run_name}_host  --out-dir ${out_dir}/host"
  fi
  log_info "  統合レポート: pic report --out-dir ${out_dir}"
}
