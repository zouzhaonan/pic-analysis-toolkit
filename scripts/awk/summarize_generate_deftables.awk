BEGIN {
  header = "count_prefix\tbarcode\tsample\tgroup"
}
{
  if (tolower($1) == "fastq_prefix" && tolower($2) == "genome" && tolower($3) == "barcode" && tolower($4) == "sample" && tolower($5) == "group") {
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
  print count_prefix, $3, $4, $5 > outfile
}
