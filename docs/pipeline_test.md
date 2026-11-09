# PIC Analysis Toolkit Test

This guide runs the basic pic workflow using a human (hg38) test dataset.
The test data are about 434,000 read pairs randomly subsampled from real PIC-seq data (4 samples in 2 groups: `R6_CTRL` and `Lig_CTRL` with *n* = 2 each); about 22 MB.
A single `pic all` command runs mapping → DESeq2 → ORA/GSEA → HTML report, and you then open the report in a browser.

## 1. Prepare the test dataset

```bash
mkdir -p "$HOME/tmp"
cd "$HOME/tmp"
curl -L "https://github.com/zouzhaonan/pic-analysis-toolkit/releases/latest/download/pic_test.tar.gz" | tar zx

wd="$HOME/tmp/pic_test"
```

The hg38 reference must be installed (see the [Installation Guide](INSTALLATION_GUIDE.md)).

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

The run is finished when `pic all: 完了しました` ("completed") is printed. It takes a few minutes (about 3 minutes on an Apple Silicon Mac) (mapping is followed by enrichment and report generation).

Check `mapping_sum`: all 4 samples should be listed, each with roughly 80,000-150,000 reads in `total`.

```bash
cat "$wd/summary/mapping_sum__test.tsv"
```

## 4. Open the report

```bash
open "$wd/report_test_hg38.html"
```

The report embeds all of its data in a single HTML file, so you can share just this file (about 6 MB at this data size).
`R6_CTRL` vs `Lig_CTRL` gives some DEGs.
The main outputs are (`<run>_<genome>` = `test_hg38`):

| Location | Contents |
| --- | --- |
| `$wd/summary/` | `mapping_sum__test.tsv`, `deftable_test_hg38.tsv`, etc. |
| `$wd/deseq2/` | `stats_test_hg38.csv`, `normalizedCountTable_test_hg38.csv`, `DEG_test_hg38/`, etc. |
| `$wd/enrich/` | `GSEA/<method>/*.csv`, `ORA/<method>/*.csv` |
| `$wd/report_test_hg38.html` | HTML report |

## 5. Advanced: run subcommands individually

To change the FDR or enrichment methods, or to run ChIP-Atlas analysis or per-gene plots, run the subcommands individually on the `pic all` results.

### 5-1. DESeq2 (changing the FDR)

```bash
pic deseq2 \
  --deftable "$wd/summary/deftable_test_hg38.tsv" \
  --count-dir "$wd/counts" \
  --genome hg38 \
  --fdr 0.05 \
  --out-dir "$wd/custom/deseq2"
```

### 5-2. ORA/GSEA (choosing methods)

```bash
pic enrich \
  --stats "$wd/custom/deseq2/stats_test_hg38.csv" \
  --genome hg38 \
  --methods GO_BP,KEGG \
  --out-dir "$wd/custom/enrich" \
  --deg-clusters "$wd/custom/deseq2/DEG_test_hg38/DEGCluster_test_hg38/DEGCluster_gene_for_ora_test_hg38.csv"
```

If `GSEA/<method>/*.csv` and `ORA/<method>/*.csv` are generated under `$wd/custom/enrich/`, the run succeeded.

### 5-3. ChIP-Atlas analysis

```bash
pic chipatlas \
  --deftable "$wd/summary/deftable_test_hg38.tsv" \
  --umi "$wd/deseq2/UMI_count_test_hg38.csv" \
  --stats "$wd/deseq2/stats_test_hg38.csv" \
  --genome hg38 \
  --out-dir "$wd/chipatlas"
```

The computation runs on DDBJ WABI (ChIP-Atlas). Results of all comparisons are collected in `<out-dir>/chipatlas_runs.csv` and `<out-dir>/chipatlas_results.csv`. Use `--comparisons` to choose comparisons and `--detach` to run in the background; rerunning with the same `--out-dir` resumes.

EA splits the DEGs of A vs B into "higher in A" and "higher in B" and uses each as the other's background (against all genes the tissue's own character dominates; against each other only features that changed between the sides remain, and features shared equally do not stand out). There is no cut-off on DEG numbers; an EA is skipped only when one side is empty. PAGE takes the raw count table and runs for every comparison. Defaults: TFs and others / Histone / ATAC-Seq, no limit on concurrent jobs.

### 5-4. plot-expression

This step plots normalized counts by group for selected genes.

```bash
cat > "$wd/genes.txt" <<'EOF'
PKM
PEG10
APOH
KRT18
EOF

pic plot-expression \
  --count "$wd/deseq2/normalizedCountTable_test_hg38.csv" \
  --deftable "$wd/summary/deftable_test_hg38.tsv" \
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
