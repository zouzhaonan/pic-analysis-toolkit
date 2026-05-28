# PIC analysis toolkit

## 概要 / Overview

PIC analysis toolkit は、PIC RNA-seq データの一次解析から統計解析・enrichment までを一連で実行できます。

PIC analysis toolkit is a command-line workflow for analyzing PIC RNA-seq data, covering primary processing through statistical and enrichment analysis.

## 注意事項 / Notes

Apple Silicon Mac でのみ動作確認しています。

Tested only on Apple Silicon Mac.

**⚠️ Intel Mac は未検証です。Windows では動きません。⚠️**

**⚠️ Intel Mac is untested. Windows is not supported. ⚠️**

## このツールでできること / What This Toolkit Can Do

- `mapping`: 一次解析を一括実行 / Run full primary workflow in one command
- `deseq2`: DESeq2 解析 / Run DESeq2 analysis
- `enrich`: ORA/GSEA 解析 / Run ORA/GSEA enrichment analysis
- `chipatlas`: ChIP-Atlas enrichment/PAGE 実行 / Run ChIP-Atlas enrichment/PAGE from DESeq2 outputs
- `plot-expression`: 指定遺伝子の発現プロット / Plot expression for selected genes

```bash
pic
pic <subcommand> --help
```

## チュートリアル / Documentations

### 1) 初回導入 / First-time setup

- 日本語: [導入方法](docs/導入方法.md)
- English: [installation_guide](docs/installation_guide.md)

### 2) テストデータで実行 / End-to-end pipeline test

- 日本語: [実行例](docs/実行例.md)
- English: [pipeline_test](docs/pipeline_test.md)

### 3) 新規ゲノム導入 / Add a new genome

- 日本語: [新規ゲノム導入](docs/新規ゲノム導入.md)
- English: [add_new_genome](docs/add_new_genome.md)
