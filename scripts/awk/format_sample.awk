#!/usr/bin/awk

# 役割:
#   整形済み sample sheet を、barcode 分割を並列化するための
#   chunk メタデータへ変換する。
# 入力:
#   stdin から受け取るタブ区切りの sample sheet。
#   Col1: FASTQ filename prefix, Col3: barcode, Col4: sample name。
#   変数: threads, raw_fastq_dir。
# 出力:
#   chunk ごとの chunk id、FASTQ prefix、read range、
#   barcode 一覧、sample 一覧。

function min(a, b) {
  return a < b ? a : b
}

function build_fastq_prefix(prefix) {
  if (prefix ~ /^\//) return prefix
  return raw_fastq_dir "/" prefix
}

NR == 1 {
  if (tolower($1) == "fastq_prefix" && tolower($2) == "genome" && tolower($3) == "barcode" && tolower($4) == "sample" && tolower($5) == "group") next
}

{
  prefix = build_fastq_prefix($1)
  if (!seen[prefix]++) {
    paths[++id] = prefix
    # 入力 FASTQ prefix ごとに一度だけリード数を数え、threads 数で分割する。
    cmd = "gzcat " prefix "_R1_001.fastq.gz | wc -l"
    while (cmd| getline lines) reads[id] = lines / 4
    chunk_size[id] = int(reads[id] / threads + 0.999999999)
  }
  barcodes[id] = barcodes[id] ? barcodes[id] "," $3 : $3
  samples[id] = samples[id] ? samples[id] "," $4 : $4
} END {
  for (i=1; i<=length(paths); i++) {
    for (c=1; c<=threads; c++) {
      start = chunk_size[i] * (c - 1) + 1
      end = min(start + chunk_size[i] - 1, reads[i])
      print c, paths[i], start ":" end, barcodes[i], samples[i]
    }
  }
}
