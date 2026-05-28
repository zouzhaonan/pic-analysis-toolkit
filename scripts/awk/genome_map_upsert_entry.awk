BEGIN { FS = "\t"; OFS = "\t" }
NR == 1 {
  print "genome", "fasta_source", "gtf_source"
  next
}
$1 == genome {
  print genome, fasta_src, gtf_src
  seen = 1
  next
}
{
  print $1, $2, $3
}
END {
  if (!seen) {
    print genome, fasta_src, gtf_src
  }
}
