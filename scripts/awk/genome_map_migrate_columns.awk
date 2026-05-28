BEGIN { FS = "\t"; OFS = "\t" }
NR == 1 {
  for (i = 1; i <= NF; i++) {
    idx[$i] = i
  }
  print "genome", "fasta_source", "gtf_source"
  next
}
{
  print \
    (idx["genome"] ? $(idx["genome"]) : ""), \
    (idx["fasta_source"] ? $(idx["fasta_source"]) : ""), \
    (idx["gtf_source"] ? $(idx["gtf_source"]) : "")
}
