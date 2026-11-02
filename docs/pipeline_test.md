# PIC Analysis Toolkit Test

This guide runs the basic pic workflow using a mouse `Brain` vs `Ovary` test dataset (*n* = 2).
A single `pic all` command runs mapping → DESeq2 → ORA/GSEA → HTML report, and you then open the report in a browser.

## 1. Prepare the test dataset

```bash
mkdir -p "$HOME/tmp"
cd "$HOME/tmp"
curl -L "https://chip-atlas.dbcls.jp/data/share/pic_test.tar.gz" | tar zx

wd="$HOME/tmp/pic_test"
```

## 2. Activate the conda environment

```bash
shell_name="$(basename "${SHELL:-zsh}")"
eval "$("$(command -v conda)" "shell.${shell_name}" hook)"
conda activate pic
```

If you forget `conda activate pic`, `pic` is either not found or, on Linux, a different program (groff's `pic`) runs.

## 3. Run everything with `pic all`

```bash
pic all \
  --sample-sheet "$wd/sample_sheet__test.txt" \
  --raw-fastq-dir "$wd/fastq" \
  --out-dir "$wd" \
  --run-name test
```

The run is finished when `pic all: 完了しました` ("completed") is printed.

Check `mapping_sum`: both `brain` and `ovary` should be listed with 2 samples each.

```bash
cat "$wd/summary/mapping_sum__test.tsv"
```

## 4. Open the report

```bash
open "$wd/report_test_mm10.html"
```

The report embeds all of its data in a single HTML file, so you can share just this file.
The main outputs are (`<run>_<genome>` = `test_mm10`):

| Location | Contents |
| --- | --- |
| `$wd/summary/` | `mapping_sum__test.tsv`, `deftable_test_mm10.tsv`, etc. |
| `$wd/deseq2/` | `stats_test_mm10.csv`, `normalizedCountTable_test_mm10.csv`, `DEG_test_mm10/`, etc. |
| `$wd/enrich/` | `GSEA/<method>/*.csv`, `ORA/<method>/*.csv` |
| `$wd/report_test_mm10.html` | HTML report |

## 5. Advanced: run subcommands individually

To change the FDR or enrichment methods, or to run ChIP-Atlas analysis or per-gene plots, run the subcommands individually on the `pic all` results.

### 5-1. DESeq2 (changing the FDR)

```bash
pic deseq2 \
  --deftable "$wd/summary/deftable_test_mm10.tsv" \
  --count-dir "$wd/counts" \
  --genome mm10 \
  --fdr 0.05 \
  --out-dir "$wd/custom/deseq2"
```

### 5-2. ORA/GSEA (choosing methods)

```bash
pic enrich \
  --stats "$wd/custom/deseq2/stats_test_mm10.csv" \
  --genome mm10 \
  --methods GO_BP,KEGG \
  --out-dir "$wd/custom/enrich" \
  --deg-clusters "$wd/custom/deseq2/DEG_test_mm10/DEGCluster_test_mm10/DEGCluster_gene_for_ora_test_mm10.csv"
```

If `GSEA/<method>/*.csv` and `ORA/<method>/*.csv` are generated under `$wd/custom/enrich/`, the run succeeded.

### 5-3. ChIP-Atlas analysis

```bash
pic chipatlas \
  --deftable "$wd/summary/deftable_test_mm10.tsv" \
  --umi "$wd/deseq2/UMI_count_test_mm10.csv" \
  --deg-list-dir "$wd/deseq2/DEG_test_mm10/DEGList_test_mm10" \
  --genome mm10 \
  --out-dir "$wd/chipatlas"
```

The computation runs on the ChIP-Atlas side. Results are saved under `<out-dir>/result/` after completion.

### 5-4. plot-expression

This step plots normalized counts by group for selected genes.

```bash
cat > "$wd/genes.txt" <<'EOF'
Map2
Gfap
Ddx4
Kit
EOF

pic plot-expression \
  --count "$wd/deseq2/normalizedCountTable_test_mm10.csv" \
  --deftable "$wd/summary/deftable_test_mm10.tsv" \
  --gene-list "$wd/genes.txt" \
  --out "$wd/plot_expression_test.pdf"
```

Warnings may appear during plotting, but this is fine as long as `$wd/plot_expression_test.pdf` is generated correctly.

## 6. Cleanup

### 6-1. Deactivate conda environment

```bash
conda deactivate
```

### 6-2. Remove test data (optional)

```bash
rm -rf "$wd"
```
