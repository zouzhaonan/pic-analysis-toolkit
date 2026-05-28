# pic Toolkit Installation Guide

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

## 2. Download the `pic` package

Unzip the downloaded `pic-channel` ZIP file, then place the extracted `pic-channel` folder anywhere you like.
The instructions below assume it is placed at `~/Downloads/pic-channel`.

**You can delete this folder after installation.**

---

## 3. Install

Run the following commands in order.

### 3-1. Remove existing `pic` environment (just in case)

If warning-like messages appear here, it is fine at this step.

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

## 4. Verify startup

### 4-1. Activate environment

```bash
shell_name="$(basename "${SHELL:-zsh}")"
eval "$("$(command -v conda)" "shell.${shell_name}" hook)"
conda activate pic
```

### 4-2. Show help

```bash
pic
```

If help is displayed, installation is complete.

## 5. Prepare references for human (hg38) and mouse (mm10)

For other species/genomes, see [Add New Genome](docs/add_new_genome.md).
Here, we place hg38 and mm10 so that analysis can start immediately.

```bash
mkdir -p $HOME/local/lib
cd $HOME/local/lib
curl https://chip-atlas.dbcls.jp/data/share/pic.tar.gz| tar zx

cd
```

## 6. Test the workflow

After completing the steps above, run the workflow check by following [Pipeline Test](pipeline_test.md).
