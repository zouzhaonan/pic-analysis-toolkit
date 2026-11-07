#!/usr/bin/env bash

# 役割:
#   `pic reanalyze`: 出来上がった解析から、サンプル (と group) を除いて deseq2 / enrich / report を
#   作り直す。mapping はやり直さない。元の出力ディレクトリは書き換えない。
# 入力:
#   --out-dir <元> (pic all / mapping の出力。summary/ の deftable・mapping_sum と counts/ bw/ ...)
#   --new-out-dir <新> (無いか空のディレクトリ)
#   --exclude-samples a,b,...  --exclude-groups g1,g2,...  (どちらか一方は必須)
#   --threads N
# 動き:
#   1. 除外を確かめる (知らないサンプル・group、group が知らないうちに 0 になる、比較が組めない は止める)。
#   2. 新しい置き場に、元の mapping の出力 (bam・bw・counts・log・demux など) への symlink を置く。
#      summary/ は写して、deftable・mapping_sum・aggregate_profile から除外したサンプルの行を外す。
#      sample_sheet.tsv も同じく外した写しにする。
#      aggregate_profile はサンプルごとの値 (bw の平均) で他のサンプルに依らないので、行を外すだけにして
#      aggregate はやり直さない。
#   3. summary/analysis_params.tsv に excluded_samples / excluded_groups / reanalyzed_from を残す
#      (元が reanalyze の出力なら、前の除外に足していく)。レポートの Overview はこれを読んで除外を示す。
#   4. pipeline_run_downstream (cli_pipeline.sh。deseq2 / enrich / report を依存関係の順に流す) を呼ぶ。
# 出力:
#   <new-out-dir>/deseq2/ enrich/ summary/ report_<run>_<genome>.html
# 注記:
#   xenograft の出力 (graft/ host/ の画分) には対応しない (止める)。

# a,b,c を 1 行に 1 つ (空は捨てる、前後の空白は取る)
reanalyze_split_csv() {
  local IFS=','
  local x
  for x in $1; do
    x="${x#"${x%%[![:space:]]*}"}"
    x="${x%"${x##*[![:space:]]}"}"
    [[ -n "$x" ]] && printf "%s\n" "$x"
  done
}

# TSV / CSV から、sample 列が除外サンプル・group 列が除外 group の行を外す。
#   reanalyze_filter_table <in> <out> <区切り (tab|comma)> <samples (改行区切り)> <groups (改行区切り)>
# sample / group の列は見出しの名前で探す (無い列は見ない)。CSV の値は引用符を外して比べる。
reanalyze_filter_table() {
  local in="$1" out="$2" sep="$3" xs="$4" xg="$5"
  local fs=$'\t'
  [[ "$sep" == comma ]] && fs=","
  awk -v FS="$fs" -v xs="$xs" -v xg="$xg" '
    function unq(v) { gsub(/\r$/, "", v); if (v ~ /^".*"$/) { v = substr(v, 2, length(v) - 2); gsub(/""/, "\"", v) } return v }
    BEGIN {
      n = split(xs, a, "\n"); for (i = 1; i <= n; i++) if (a[i] != "") S[a[i]] = 1
      m = split(xg, b, "\n"); for (i = 1; i <= m; i++) if (b[i] != "") G[b[i]] = 1
    }
    NR == 1 { for (i = 1; i <= NF; i++) { h = unq($i); if (h == "sample") sc = i; if (h == "group") gc = i } print; next }
    {
      if (sc && (unq($sc) in S)) next
      if (gc && (unq($gc) in G)) next
      print
    }' "$in" > "$out"
}

# 除外を deftable に当てて確かめる。問題があればメッセージを出して 1 を返す。
#   reanalyze_check_deftable <deftable> <samples> <groups>
reanalyze_check_deftable() {
  local f="$1" xs="$2" xg="$3"
  awk -F'\t' -v xs="$xs" -v xg="$xg" -v f="$(basename "$f")" '
    BEGIN {
      n = split(xs, a, "\n"); for (i = 1; i <= n; i++) if (a[i] != "") S[a[i]] = 1
      m = split(xg, b, "\n"); for (i = 1; i <= m; i++) if (b[i] != "") G[b[i]] = 1
    }
    NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "sample") sc = i; if ($i == "group") gc = i } next }
    {
      gsub(/\r$/, "")
      s = $sc; g = $gc
      if (!(g in tot)) order[++ng] = g
      tot[g]++
      if (!((s in S) || (g in G))) { keep[g]++; nkeep++ }
    }
    END {
      if (!sc || !gc) { printf "[ERROR] %s: sample / group 列がありません\n", f; exit 1 }
      bad = 0; left = 0
      for (k = 1; k <= ng; k++) {
        g = order[k]
        if (keep[g] > 0) { left++; continue }
        if (!(g in G)) {
          printf "[ERROR] %s: group \"%s\" のサンプルがすべて除外されます (0 になる)。group ごと外すなら --exclude-groups %s を付けてください\n", f, g, g
          bad = 1
        }
      }
      if (left < 2) {
        printf "[ERROR] %s: 除外の後に残る group が %d 個です。比較 (2 group 以上) が組めません\n", f, left
        bad = 1
      }
      if (!bad) printf "[INFO] %s: 残り %d サンプル / %d group\n", f, nkeep, left
      exit bad
    }' "$f"
}

# パスを、作らずに絶対パス (symlink を解いたもの) にする。存在する一番深い親を pwd -P で解き、残りを足す。
reanalyze_abs_path() {
  local p="$1" rest="" b
  [[ "$p" == /* ]] || p="$(pwd -P)/${p}"
  while [[ ! -d "$p" ]]; do
    b="$(basename "$p")"
    rest="/${b}${rest}"
    p="$(dirname "$p")"
  done
  printf "%s%s\n" "$(cd "$p" && pwd -P)" "$rest"
}

run_reanalyze_subcommand() {
  if [[ $# -eq 0 ]]; then
    show_help_for_subcommand reanalyze
    return 0
  fi
  local src_in="" dst_in="" xs_csv="" xg_csv="" threads="$PIC_DEFAULT_THREADS"
  while [[ $# -gt 0 ]]; do
    case "$1" in
    --out-dir) src_in="${2:-}"; shift 2 ;;
    --new-out-dir) dst_in="${2:-}"; shift 2 ;;
    --exclude-samples) xs_csv="${2:-}"; shift 2 ;;
    --exclude-groups) xg_csv="${2:-}"; shift 2 ;;
    --threads) threads="${2:-}"; shift 2 ;;
    --help|-h) show_help_for_subcommand reanalyze; return 0 ;;
    *)
      handle_error "Unknown option: $1"
      handle_error "pic reanalyze のオプション: --out-dir / --new-out-dir / --exclude-samples / --exclude-groups / --threads"
      exit 1
      ;;
    esac
  done

  [[ -n "$src_in" ]] || { handle_error "--out-dir (元の出力ディレクトリ) を指定してください"; exit 1; }
  [[ -n "$dst_in" ]] || { handle_error "--new-out-dir (新しい出力ディレクトリ) を指定してください"; exit 1; }
  if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
    handle_error "--threads は 1 以上の整数で指定してください: ${threads}"
    exit 1
  fi
  local xs xg
  xs="$(reanalyze_split_csv "$xs_csv" | awk '!seen[$0]++')"
  xg="$(reanalyze_split_csv "$xg_csv" | awk '!seen[$0]++')"
  if [[ -z "$xs" && -z "$xg" ]]; then
    handle_error "--exclude-samples か --exclude-groups のどちらかで、除外するものを指定してください"
    exit 1
  fi

  # --- 元の置き場 ---
  [[ -d "$src_in" ]] || { handle_error "--out-dir がありません: ${src_in}"; exit 1; }
  local src
  src="$(cd "$src_in" && pwd -P)"
  if [[ ! -d "${src}/summary" ]]; then
    handle_error "${src}/summary がありません (pic all / mapping の出力ディレクトリを指定してください)"
    exit 1
  fi
  if compgen -G "${src}/xenograft_classify_summary__*.tsv" >/dev/null || [[ -d "${src}/graft" || -d "${src}/host" ]]; then
    handle_error "xenograft の出力には pic reanalyze は対応していません: ${src}"
    exit 1
  fi
  local msums=() f
  for f in "${src}"/summary/mapping_sum__*.tsv; do [[ -e "$f" ]] && msums+=("$f"); done
  if [[ "${#msums[@]}" -ne 1 ]]; then
    handle_error "${src}/summary に mapping_sum__<run>.tsv が 1 つだけあるはずですが ${#msums[@]} 個です"
    exit 1
  fi
  local run
  run="$(basename "${msums[0]}" .tsv)"
  run="${run#mapping_sum__}"
  local deftables=() genomes=() g
  for f in "${src}"/summary/deftable_"${run}"_*.tsv; do
    [[ -e "$f" ]] || continue
    deftables+=("$f")
    g="$(basename "$f" .tsv)"
    genomes+=("${g#deftable_"${run}"_}")
  done
  if [[ "${#deftables[@]}" -eq 0 ]]; then
    handle_error "deftable が見つかりません: ${src}/summary/deftable_${run}_*.tsv"
    exit 1
  fi

  # --- 除外を確かめる (何も作る前に) ---
  local known_s known_g x bad=0
  known_s="$(for f in "${deftables[@]}"; do awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="sample") c=i; next} c{gsub(/\r$/,""); print $c}' "$f"; done | sort -u)"
  known_g="$(for f in "${deftables[@]}"; do awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="group") c=i; next} c{gsub(/\r$/,""); print $c}' "$f"; done | sort -u)"
  while IFS= read -r x; do
    [[ -z "$x" ]] && continue
    if ! grep -Fxq -- "$x" <<<"$known_s"; then handle_error "知らないサンプルです: ${x} (deftable にありません)"; bad=1; fi
  done <<<"$xs"
  while IFS= read -r x; do
    [[ -z "$x" ]] && continue
    if ! grep -Fxq -- "$x" <<<"$known_g"; then handle_error "知らない group です: ${x} (deftable にありません)"; bad=1; fi
  done <<<"$xg"
  [[ "$bad" = 0 ]] || exit 1
  for f in "${deftables[@]}"; do
    reanalyze_check_deftable "$f" "$xs" "$xg" >&2 || bad=1
  done
  [[ "$bad" = 0 ]] || { handle_error "除外を当てると解析が組めません。何も作らずに止めます"; exit 1; }

  # --- 新しい置き場 ---
  if [[ -e "$dst_in" && ! -d "$dst_in" ]]; then
    handle_error "--new-out-dir がディレクトリではありません: ${dst_in}"
    exit 1
  fi
  if [[ -d "$dst_in" ]] && [[ -n "$(ls -A "$dst_in" 2>/dev/null)" ]]; then
    handle_error "--new-out-dir が空ではありません (上書きしません): ${dst_in}"
    exit 1
  fi
  local dst
  dst="$(reanalyze_abs_path "$dst_in")"
  if [[ "$dst" == "$src" || "$dst" == "$src"/* || "$dst" == */.. || "$dst" == */../* || "$dst" == */. || "$dst" == */./* ]]; then
    handle_error "--new-out-dir は元の --out-dir の外に、. や .. を含まない形で指定してください (元を書き換えないため): ${dst}"
    exit 1
  fi
  mkdir -p "$dst" || { handle_error "--new-out-dir を作れません: ${dst}"; exit 1; }
  log_info "pic reanalyze: ${src} -> ${dst} (run=${run}, genome=${genomes[*]})"

  # mapping の出力は symlink。作り直すもの (summary・deseq2・enrich・report・tmp) と sample_sheet は除く。
  local e b
  for e in "${src}"/* "${src}"/.[!.]*; do
    [[ -e "$e" || -L "$e" ]] || continue
    b="$(basename "$e")"
    case "$b" in
      summary|deseq2|enrich|tmp|report_*.html|sample_sheet*.tsv) continue ;;
    esac
    ln -s "${src}/${b}" "${dst}/${b}" || { handle_error "symlink を作れません: ${dst}/${b}"; exit 1; }
  done

  # summary/ は実体を写して (元への symlink にしない)、サンプルごとの表から除外した行を外す。
  mkdir -p "${dst}/summary"
  for f in "${src}"/summary/*; do
    [[ -f "$f" ]] || continue
    b="$(basename "$f")"
    case "$b" in
      deftable_*.tsv|mapping_sum__*.tsv) reanalyze_filter_table "$f" "${dst}/summary/${b}" tab "$xs" "$xg" ;;
      aggregate_profile_*.csv) reanalyze_filter_table "$f" "${dst}/summary/${b}" comma "$xs" "$xg" ;;
      *) cp -p "$f" "${dst}/summary/${b}" ;;
    esac
  done
  for f in "${src}"/sample_sheet*.tsv; do
    [[ -f "$f" ]] || continue
    reanalyze_filter_table "$f" "${dst}/$(basename "$f")" tab "$xs" "$xg"
  done

  # 除外の記録 (前に reanalyze した出力からなら、前の除外に足す)
  local params="${dst}/summary/analysis_params.tsv" prev_s="" prev_g="" all_s all_g
  if [[ -f "$params" ]]; then
    prev_s="$(awk -F'\t' '$1=="excluded_samples"{print $2}' "$params")"
    prev_g="$(awk -F'\t' '$1=="excluded_groups"{print $2}' "$params")"
    awk -F'\t' '$1!="excluded_samples" && $1!="excluded_groups" && $1!="reanalyzed_from"' "$params" > "${params}.tmp"
    mv "${params}.tmp" "$params"
  else
    printf "key\tvalue\n" > "$params"
  fi
  all_s="$( { reanalyze_split_csv "$prev_s"; printf "%s\n" "$xs"; } | awk 'NF && !seen[$0]++' | paste -sd, -)"
  all_g="$( { reanalyze_split_csv "$prev_g"; printf "%s\n" "$xg"; } | awk 'NF && !seen[$0]++' | paste -sd, -)"
  {
    [[ -n "$all_s" ]] && printf "excluded_samples\t%s\n" "$all_s"
    [[ -n "$all_g" ]] && printf "excluded_groups\t%s\n" "$all_g"
    printf "reanalyzed_from\t%s\n" "$src"
  } >> "$params"

  # --- deseq2 / enrich / report (aggregate は行を外した写しを使うのでやり直さない) ---
  declare -g OUTPUT_DIR="$dst"
  declare -g RUN_NAME="$run"
  declare -g DEFTABLE_DIR="${dst}/summary"
  declare -g COUNTS_DIR="${dst}/counts"
  declare -g THREADS="$threads"
  pipeline_run_downstream 1 0 "${genomes[@]}"

  local missing=0
  for g in "${genomes[@]}"; do
    if [[ ! -s "${dst}/report_${run}_${g}.html" ]]; then
      handle_error "pic reanalyze: レポートができていません: ${dst}/report_${run}_${g}.html"
      missing=1
    fi
  done
  [[ "$missing" = 0 ]] || exit 1
  log_info "pic reanalyze: 完了しました (new-out-dir=${dst})"
}
