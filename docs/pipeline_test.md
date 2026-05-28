# PIC Analysis Toolkit Test

This guide runs the basic pic workflow using a mouse `Brain` vs `Ovary` test dataset (*n* = 2).

## 1. Prepare the test dataset

```bash
mkdir -p "$HOME/tmp"
cd "$HOME/tmp"
curl -L "https://chip-atlas.dbcls.jp/data/share/pic_test.tar.gz" | tar zx

wd="$HOME/tmp/pic_test"
```

## 2. Run mapping

### 2-1. Activate the conda environment

```bash
shell_name="$(basename "${SHELL:-zsh}")"
eval "$("$(command -v conda)" "shell.${shell_name}" hook)"
conda activate pic
```

### 2-2. Run mapping

```bash
pic mapping \
  --sample-sheet "$wd/sample_sheet__test.txt" \
  --raw-fastq-dir "$wd/fastq" \
  --out-dir "$wd" \
  --run-name test
```

After completion, check `mapping_sum`:

```bash
cat "$wd/mapping_sum__test.tsv"
```

If both `brain` and `ovary` are shown with 2 samples each, mapping succeeded.

## 3. Run DESeq2

`deftable` is automatically generated after `pic mapping`.

```bash
pic deseq2 \
  --deftable "$wd/deftable_mm10_test.tsv" \
  --count-dir "$wd/counts" \
  --genome mm10 \
  --out-dir "$wd/deseq2"
```

After completion, check `stats`:

```bash
cat "$wd/deseq2/stats_mm10_test.csv"
```

## 4. Run ORA/GSEA enrichment

```bash
pic enrich \
  --stats "$wd/deseq2/stats_mm10_test.csv" \
  --genome mm10 \
  --out-dir "$wd/enrich" \
  --deg-clusters "$wd/deseq2/DEG/DEGCluster/DEGCluster_gene_for_ora_mm10_test.csv"
```

If `plots/` and `csv/` are generated under `$wd/enrich`, the run succeeded.

## 5. Run ChIP-Atlas analysis

```bash
pic chipatlas \
  --deftable "$wd/deftable_mm10_test.tsv" \
  --umi "$wd/deseq2/UMI_count_mm10_test.csv" \
  --deg-list-dir "$wd/deseq2/DEG/DEGList" \
  --genome mm10 \
  --out-dir "$wd/chipatlas"
```

The computation runs on the ChIP-Atlas side. Results are saved under the specified `out-dir` after completion.

## 6. Run plot-expression

This step plots normalized counts by group for selected genes.

First, create a gene list for plotting.

```bash
cat > "$wd/genes.txt" <<'EOF'
Map2
Gfap
Ddx4
Kit
EOF
```

Then generate expression plots from the DESeq2 normalized count table.

```bash
pic plot-expression \
  --count "$wd/deseq2/normalizedCountTable_mm10_test.csv" \
  --deftable "$wd/deftable_mm10_test.tsv" \
  --gene-list "$wd/genes.txt" \
  --out "$wd/deseq2/plot_expression_test.pdf"
```

Warnings may appear during plotting, but this is fine as long as `$wd/deseq2/plot_expression_test.pdf` is generated correctly.

## 7. Cleanup

### 7-1. Deactivate conda environment

```bash
conda deactivate
```

### 7-2. Remove test data (optional)

```bash
rm -rf "$wd"
```
