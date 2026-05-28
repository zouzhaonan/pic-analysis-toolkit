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
  print $0, "all"
  for (point in points) {
    if (nreads[$4] / point >= 3) print $0, point
  }
}
