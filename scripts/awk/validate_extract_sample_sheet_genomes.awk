NR == 1 {
  first = tolower($1)
  second = tolower($2)
  if (first == "fastq_prefix" && second == "genome") next
}
$2 != "" {
  gsub(/\r/, "", $2)
  print $2
}
