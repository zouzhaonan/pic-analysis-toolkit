NR == 1 {
  gsub(/\\t/, "\t", $0)
  print
  next
}
{ print }
