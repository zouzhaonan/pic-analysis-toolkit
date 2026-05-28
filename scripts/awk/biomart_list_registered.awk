BEGIN { FS = "\t" }
NR == 1 {
  for (i = 1; i <= NF; i++) {
    h = $i
    gsub(/\r/, "", h)
    if (h == "genome") g = i
    if (h == "biomart_dataset") d = i
  }
  next
}
g > 0 && d > 0 {
  gv = $g
  dv = $d
  gsub(/\r/, "", gv)
  gsub(/\r/, "", dv)
  if (gv != "" && dv != "" && gv != "NA" && dv != "NA") {
    print gv "\t" dv
  }
}
