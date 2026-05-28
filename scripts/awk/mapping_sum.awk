#!/usr/bin/awk

# 役割:
#   mapping log、featureCounts の結果、UMI count をまとめて、
#   sample ごとの集計結果を1行ずつ作る。
# 入力:
#   stdin から受け取る total_reads table。
#   変数: job_queue, tmp_dir, count_dir。
# 出力:
#   mapping と UMI の指標をまとめたタブ区切りの集計行。

BEGIN {
  while ((getline < job_queue) > 0) {
    barcode[$4] = $3
    group[$4] = $5
    points[$4] = points[$4] ? points[$4] "," $6 : $6
  }
  close(job_queue)
} {
  sample = $1
  split(points[sample], point_arr, ",")

  for (i=1; i<=length(point_arr); i++) {
    point = point_arr[i]
    total = point == "all" ? $2 : point

    for (key in reads) delete reads[key]
    trimmed = multi = 0
    umis = genes = 0

    # simulation では sample_point、通常モードでは sample を使う。
    stem = (point == "all") ? sample : sample "_" point
    mapping_log = tmp_dir "/" stem ".mapping.log"
    feature_sum = tmp_dir "/" stem ".feature.tsv.summary"
    umi_count = count_dir "/" stem ".txt"
    umis = genes = 0

    while ((getline < mapping_log) > 0) {
      if ($0 ~ " reads; of these:") {
        gsub(" reads; of these:", "")
        trimmed = total - $0
      } else if ($0 ~ "aligned >1 times") {
        split($0, parts, " \\(")
        gsub(/ +/, "", parts[1])
        multi = parts[1]
      }
    }

    while ((getline < feature_sum) > 0) {
      gsub("Unassigned_", "", $1)
      reads[$1] = $2
    }

    while ((getline < umi_count) > 0) {
      if ($3 > 0) {
        umis += $3
        genes++
      }
    }

    close(mapping_log)
    close(feature_sum)
    close(umi_count)
 
    # 追加の指標として、1遺伝子あたりの UMI 数などを計算する。
    umi_per_gene = genes > 0 ? umis/genes : "NA"
    ass_per_umis = umis > 0 ? reads["Assigned"]/umis : "NA"

    print sample, group[sample], barcode[sample], total, trimmed, reads["Unmapped"], multi, \
          reads["NoFeatures"], reads["Ambiguity"], reads["Assigned"], \
          umis, genes, umi_per_gene, ass_per_umis
  }

}
