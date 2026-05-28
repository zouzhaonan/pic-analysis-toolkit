#!/usr/bin/awk

# 役割:
#   連結済みの R1/R2 FASTQ を barcode ごとに振り分け、
#   sample ごとの FASTQ 断片として書き出す。
#   あわせて read 名に UMI 情報を残す。
# 入力:
#   stdin から受け取る連結済み paired FASTQ。
#   変数: barcodes, samples, chunk, tmp_dir。
# 出力:
#   tmp_dir 配下の sample ごとの FASTQ 断片と read 数ファイル。

BEGIN {
  split(barcodes, barcode_list, ",")
  split(samples, sample_list, ",")

  for (i=1; i<=length(barcode_list); i++) {
    barcode_to_sample[barcode_list[i]] = sample_list[i]
    sample_seen[sample_list[i]] = 1
  }
}

function file_name(sample, suffix) {
  return sprintf("%s/%s.%s.%s", tmp_dir, sample, suffix, chunk)
} {
  if (NR % 4 == 1) {
    read1_header = $3
    read2_header = $4
    output_file = ""
    emit_read = 0
  } else if (NR % 4 == 2) {
    # PIC の read1 先頭 6 塩基が UMI、続く 6 塩基が barcode という前提。
    barcode = substr($1, 7, 6)
    if (barcode in barcode_to_sample) {
      umi = substr($1, 1, 6)
      sample = barcode_to_sample[barcode]
      output_file = file_name(sample, "fastq")
      nreads[sample]++
      emit_read = 1
      print read1_header "_" barcode "_" umi, read2_header > output_file
    }
  }
  if (emit_read == 1) print $2 > output_file
} END {
  # 後で chunk ごとの read 数を合計できるよう、sample ごとに保存する。
  for (sample in sample_seen) print nreads[sample] > file_name(sample, "nread")
}
