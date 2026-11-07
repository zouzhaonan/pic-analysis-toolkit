# PIC analysis toolkit

## 概要 / Overview

PIC analysis toolkit は、PIC RNA-seq データの一次解析から統計解析・enrichment・HTML レポートまでを一連で実行できます。

PIC analysis toolkit is a command-line workflow for analyzing PIC RNA-seq data, covering primary processing through statistical and enrichment analysis and an HTML report.

```bash
pic all --sample-sheet sample.tsv --raw-fastq-dir fastq/ --out-dir out/ --run-name test
```

## 注意事項 / Notes

- Apple Silicon Mac (osx-arm64) と Linux (x86_64) で動作確認しています。/ Tested on Apple Silicon Mac (osx-arm64) and Linux (x86_64).
- **⚠️ 1.5.0 から OS ごとのパッケージになり、Intel Mac 用はありません。Windows では動きません。/ Since 1.5.0 the package is built per OS; there is no Intel Mac build. Windows is not supported. ⚠️**
- conda が必要です（[Miniforge](https://github.com/conda-forge/miniforge) を推奨）。/ Requires conda ([Miniforge](https://github.com/conda-forge/miniforge) recommended).
- ネット接続が必要です（インストール、リファレンスの取得、`pic chipatlas`、`pic enrich --prepare-libs` などで外部のサーバーにつなぎます）。/ Requires an internet connection (installation, reference download, `pic chipatlas`, `pic enrich --prepare-libs`, etc. access external servers).
- `conda activate pic` を忘れると、Linux では groff の `pic` が動きます。/ Without `conda activate pic`, Linux runs groff's `pic` instead.

## ダウンロードとインストール / Download and Install

[Download Page](https://github.com/zouzhaonan/pic-analysis-toolkit/releases/latest) から最新版の `pic-channel.zip` をダウンロードし、解凍後に [導入方法](docs/導入方法.md) にしたがってインストールしてください。

Download the latest `pic-channel.zip` from the [Download Page](https://github.com/zouzhaonan/pic-analysis-toolkit/releases/latest), unzip it, then follow the [Installation Guide](docs/INSTALLATION_GUIDE.md).

## 出力 / Output

`pic all` は `--out-dir` に次を作ります（`<run>_<genome>` はゲノムごと）。/ `pic all` writes the following to `--out-dir` (one `<run>_<genome>` per genome).

| 場所 / Location | 中身 / Contents |
| --- | --- |
| `summary/` | mapping_sum, deftable, aggregate_profile |
| `deseq2/` | DESeq2 の結果 / DESeq2 results (`stats_<run>_<genome>.csv`, `DEG_<run>_<genome>/` ...) |
| `enrich/` | ORA/GSEA の結果 / ORA/GSEA results (`ORA/<method>/`, `GSEA/<method>/`) |
| `report_<run>_<genome>.html` | 自己完結の HTML レポート / self-contained HTML report |

レポートの RIN タブは、TapeStation のデータを `pic report --rin-dir` か環境変数 `PIC_RIN_DIR` で渡したときだけ出ます。/ The RIN tab appears only when TapeStation data is given via `pic report --rin-dir` or the `PIC_RIN_DIR` environment variable.

## チュートリアル / Documentation

### 1) セットアップ / Setup

- 日本語: [導入方法](docs/導入方法.md)
- English: [Installation Guide](docs/INSTALLATION_GUIDE.md)

### 2) テストデータで実行 / End-to-end pipeline test

`pic all` 1 本で最後まで実行し、HTML レポートを開きます。/ Run everything with a single `pic all` and open the HTML report.

- 日本語: [実行例](docs/実行例.md)
- English: [Pipeline Test](docs/pipeline_test.md)

### 3) 新規ゲノム導入 / Add a new genome

- 日本語: [新規ゲノム導入](docs/新規ゲノム導入.md)
- English: [Add New Genome](docs/ADD_NEW_GENOME.md)
