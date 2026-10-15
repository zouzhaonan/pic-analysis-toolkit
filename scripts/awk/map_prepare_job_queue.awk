BEGIN {
  n_points = 0
  if (sim_mode == 1) {
    while ((getline line < total_reads) > 0) {
      split(line, row, "\t")
      nreads[row[1]] = row[2]
    }
    close(total_reads)
    n_points = split("10000,100000,1000000,10000000,20000,200000,2000000,20000000,50000,500000,5000000,50000000", tmp, ",")
    for (i = 1; i <= n_points; i++) points[tmp[i]]++
  }
}
{
  if (NR == 1 && tolower($1) == "fastq_prefix" && tolower($2) == "genome" && tolower($3) == "barcode" && tolower($4) == "sample" && tolower($5) == "group") next
  # 同じ sample が複数行 (= 複数 barcode を 1 サンプルに合算) でも、demux は 1 つの
  # FASTQ にまとめるので mapping は 1 回でよい。barcode 列はカンマ連結で持ち回す。
  if ($4 in row) {
    if (("," bcs[$4] ",") !~ ("," $3 ",")) bcs[$4] = bcs[$4] "," $3
    next
  }
  row[$4] = $0
  bcs[$4] = $3
  order[++n] = $4
}
END {
  for (i = 1; i <= n; i++) {
    s = order[i]
    nf = split(row[s], f, "\t")
    f[3] = bcs[s]
    line = f[1]
    for (j = 2; j <= nf; j++) line = line OFS f[j]
    print line, "all"
    for (point in points) {
      if (nreads[s] / point >= 3) print line, point
    }
  }
}
