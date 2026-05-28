# pic Toolkit Installation Guide

Tested only on Apple Silicon Mac.
**⚠️ Intel Mac is untested. Windows is not supported. ⚠️**

If you complete this guide, you will install the toolkit and set up hg38/mm10 references.
**⚠️ The full process can take a few hours. Please run it when your network is stable and you have enough time. ⚠️**

## 1. Prerequisites

### 1-1. Install Miniconda (only if not installed)

Skip this step if `conda` is already available.

- [Download the installer](https://www.anaconda.com/download/success)
- Install it like a normal application.

### 1-2. Check that `conda` works in Terminal

```bash
conda --version
```

---

## 2. Download `pic-channel`

Place the shared `pic-channel` folder anywhere you like.
Example: `~/Downloads/pic-channel`

**You can delete this folder after installation.**

---

## 3. Install

Run the following commands in order.

### 3-1. Remove existing `pic` environment (just in case)

If you see warning-like messages here, it is fine at this step.

```bash
conda deactivate || true
conda env remove -n pic -y || true
```

### 3-2. Get the `pic-channel` path

```bash
cd ~/Downloads/pic-channel
channel_dir="$(pwd)"
```

### 3-3. Run installation

The solver may appear to be stuck for tens of minutes.
This is normal, so please wait.

```bash
conda config --set channel_priority strict

conda create -n pic -y \
  --override-channels \
  -c "file://$channel_dir" \
  -c conda-forge -c bioconda \
  pic
```

---

## 4. Verify

### 4-1. Activate environment

```bash
conda activate pic
```

### 4-2. Show help

```bash
pic
```

If help is displayed, installation is complete.

## 5. Prepare References for Human (hg38) and Mouse (mm10)

For other species/genomes, see the operation guide.
Here, we prepare hg38 and mm10 so you can start running analyses.

**⚠️ This section takes time. Run it with a stable internet connection. ⚠️**

First, move to a working directory:

```bash
tmp_dir="$HOME/Downloads/pic_tmp"
mkdir -p "$tmp_dir"
cd "$tmp_dir"
```

### 5-1. build-genome

This step prepares HISAT2 index and GTF.
Copy and run the block below as-is:

```bash
prepare_genome() {
  genome="$1"
  fasta_url="$2"
  gtf_url="$3"

  fasta_gz="$(basename "$fasta_url")"
  gtf_gz="$(basename "$gtf_url")"
  fasta="${fasta_gz%.gz}"
  gtf="${gtf_gz%.gz}"

  curl -L "$fasta_url" -o "$fasta_gz"
  curl -L "$gtf_url" -o "$gtf_gz"
  gunzip -f "$fasta_gz" "$gtf_gz"

  pic build-genome \
    --genome "$genome" \
    --fasta "$fasta" \
    --gtf "$gtf" \
    --fasta-source "$fasta_url" \
    --gtf-source "$gtf_url"
}

prepare_genome hg38 \
  "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
  "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"

prepare_genome mm10 \
  "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz" \
  "https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
```

### 5-2. Register biomart

This adds annotation used in downstream count/stat tables.

```bash
pic manage-biomart --register --genome hg38 --dataset hsapiens_gene_ensembl
pic manage-biomart --register --genome mm10 --dataset mmusculus_gene_ensembl
```

### 5-3. Prepare enrichment libraries

```bash
pic enrich --prepare-libs --genome hg38
pic enrich --prepare-libs --genome mm10
```

Notes:

- `manage-biomart` and `enrich --prepare-libs` access external servers.
- They may fail when the server is busy or temporarily unavailable.
- If failed, rerun only the relevant command(s) in [5-2](#5-2-register-biomart) or [5-3](#5-3-prepare-enrichment-libraries).

### 5-4. Remove temporary files

```bash
cd
rm -rf "$tmp_dir"
```
