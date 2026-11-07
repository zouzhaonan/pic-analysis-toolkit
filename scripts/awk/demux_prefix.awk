#!/usr/bin/awk

# 役割:
#   1 つの FASTQ の組 (R1/R2 を paste でつないだもの) を 1 回だけ読み、barcode ごとに
#   サンプルへ振り分けて、サンプルごとの pigz へ直接書く。read 名に barcode と UMI を残す。
# 入力:
#   stdin: paste <(R1) <(R2)。1 行目は "@R1名 R1注記<TAB>@R2名 R2注記"。
#   変数: barcodes, samples (カンマ区切り、同じ順), part_suffix (例 ".part.1_S1.fastq.gz"),
#         out_dir, nread_file, zthreads, barcode_first。
# 出力:
#   <out_dir>/<sample><part_suffix>.tmp (gzip。読みが 1 つも無いサンプルは作らない)、
#   nread_file (sample<TAB>読みの数。sheet にある全サンプル、0 も書く)。
#   振り分けの規則と出力の中身は、旧 split_to_chunks.awk (seqkit range で 8 分割) と同じ。

function shq(s) {
  gsub(/'/, "'\\''", s)
  return "'" s "'"
}

BEGIN {
  nb = split(barcodes, barcode_list, ",")
  split(samples, sample_list, ",")
  for (i = 1; i <= nb; i++) {
    barcode_to_sample[barcode_list[i]] = sample_list[i]
    if (!(sample_list[i] in nreads)) {
      nreads[sample_list[i]] = 0
      order[++ns] = sample_list[i]
    }
  }
}

{
  if (NR % 4 == 1) {
    read1_header = $3
    read2_header = $4
    emit_read = 0
  } else if (NR % 4 == 2) {
    # 既定は PIC 標準の read1 = UMI(6) + barcode(6)。barcode_first=1 は逆の並び。
    if (barcode_first == 1) {
      barcode = substr($1, 1, 6)
      umi = substr($1, 7, 6)
    } else {
      barcode = substr($1, 7, 6)
      umi = substr($1, 1, 6)
    }
    if (barcode in barcode_to_sample) {
      sample = barcode_to_sample[barcode]
      if (!(sample in cmd)) {
        cmd[sample] = "pigz -p " zthreads " > " shq(out_dir "/" sample part_suffix ".tmp")
      }
      out = cmd[sample]
      nreads[sample]++
      emit_read = 1
      print read1_header "_" barcode "_" umi, read2_header | out
    }
  }
  if (emit_read == 1) print $2 | out
}

END {
  # 4 行ごとに 1 リード。割り切れないのは入力が途中で切れたとき。
  if (NR % 4 != 0) {
    printf("[ERROR] truncated FASTQ input (%d lines)\n", NR) > "/dev/stderr"
    exit 2
  }
  for (s in cmd) {
    if (close(cmd[s]) != 0) {
      printf("[ERROR] pigz failed for %s\n", s) > "/dev/stderr"
      exit 2
    }
  }
  for (i = 1; i <= ns; i++) printf("%s\t%d\n", order[i], nreads[order[i]]) > nread_file
}
