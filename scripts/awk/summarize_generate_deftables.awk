# サンプルシート -> genome ごとの deftable (count_prefix, barcode, sample, group [, 追加の列…])。
# 見出し行 (fastq_prefix genome barcode sample group …) があり、group の後ろに見出し付きの列
# (例: donor) があれば、その列を deftable にも同じ見出しで写す (pic deseq2 --covariate で使う)。
# 見出しの無い列・見出し行の無いシートの 6 列目以降は今まで通り写さない。
BEGIN {
  header = "count_prefix\tbarcode\tsample\tgroup"
  n_extra = 0
}
{
  if (tolower($1) == "fastq_prefix" && tolower($2) == "genome" && tolower($3) == "barcode" && tolower($4) == "sample" && tolower($5) == "group") {
    if (!header_seen) {
      header_seen = 1
      for (i = 6; i <= NF; i++) {
        name = $i
        gsub(/^[ \t]+|[ \t]+$/, "", name)
        if (name == "") continue
        n_extra++
        extra_col[n_extra] = i
        header = header "\t" name
      }
    }
    next
  }
  genome = $2
  gsub(/[^[:alnum:]_.-]+/, "_", genome)
  outfile = deftable_dir "/deftable_" run_name "_" genome ".tsv"
  if (!(outfile in header_written)) {
    print header > outfile
    header_written[outfile] = 1
  }
  if (sim_mode == 1) {
    count_prefix = $4 "_all"
  } else {
    count_prefix = $4
  }
  line = count_prefix OFS $3 OFS $4 OFS $5
  for (k = 1; k <= n_extra; k++) {
    v = (extra_col[k] <= NF) ? $(extra_col[k]) : ""
    gsub(/^[ \t]+|[ \t]+$/, "", v)
    line = line OFS v
  }
  print line > outfile
}
