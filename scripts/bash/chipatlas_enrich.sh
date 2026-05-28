curl -X POST \
  -d "format=text" \
  -d "result=www" \
  -d "genome=mm10" \
  -d "antigenClass=TFs and others" \
  -d "cellClass=All cell types" \
  -d "threshold=50" \
  -d "typeA=gene" \
  --data-urlencode "bedAFile@brain_vs_ovary__brain.txt" \
  -d "typeB=gene" \
  --data-urlencode "bedBFile@brain_vs_ovary__ovary.txt" \
  -d "permTime=1" \
  -d "title=brain_vs_ovary" \
  -d "descriptionA=brain" \
  -d "descriptionB=ovary" \
  -d "distanceUp=5000" \
  -d "distanceDown=5000" \
  -d "sbatchOptions=-p epyc -t 180" \
  https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/

while :; do
  status="$(curl -s https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/wabi_chipatlas_2026-0522-1315-30-963-092669 | awk '$1 == "status:" {print $2}')"
  if [[ $status == "finished" ]]; then
    curl -s "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/wabi_chipatlas_2026-0522-1315-30-963-092669?info=result&format=tsv" >brain_vs_ovary__chipatlas_ea.tsv
    curl -s "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/wabi_chipatlas_2026-0522-1315-30-963-092669?info=result&format=html" >brain_vs_ovary__chipatlas_ea.html
    curl -s "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/wabi_chipatlas_2026-0522-1315-30-963-092669?info=result&format=log" >brain_vs_ovary__chipatlas_ea.log
    break
  fi
  sleep 60
done

cut -d',' -f2,5- UMI_count_mm10_qH0MSu2b.csv >UMI_count_mm10_for_chipatlas.csv

curl -X POST \
  -d "format=text" \
  -d "result=www" \
  -d "genome=mm10" \
  -d "antigenClass=TFs and others" \
  -d "cellClass=All cell types" \
  -d "threshold=50" \
  -d "typeA=count" \
  --data-urlencode "bedAFile@UMI_count_mm10_for_chipatlas.csv" \
  -d "typeB=empty" \
  -d "bedBFile=empty" \
  -d "permTime=1" \
  -d "title=brain_vs_ovary" \
  -d "descriptionA=empty" \
  -d "descriptionB=empty" \
  -d "distanceUp=5000" \
  -d "distanceDown=5000" \
  -d "sbatchOptions=-p epyc -t 180" \
  https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/

while :; do
  status="$(curl -s https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/wabi_chipatlas_2026-0522-1333-39-460-954573 | awk '$1 == "status:" {print $2}')"
  if [[ $status == "finished" ]]; then
    curl -s "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/wabi_chipatlas_2026-0522-1333-39-460-954573?info=result&format=tsv" >brain_vs_ovary__chipatlas_page.tsv
    curl -s "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/wabi_chipatlas_2026-0522-1333-39-460-954573?info=result&format=html" >brain_vs_ovary__chipatlas_page.html
    curl -s "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/wabi_chipatlas_2026-0522-1333-39-460-954573?info=result&format=log" >brain_vs_ovary__chipatlas_page.log
    break
  fi
  sleep 60
done
