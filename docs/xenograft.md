# xenograft (xengsort) 解析ガイド

xenograft サンプル（host と graft の 2 種が 1 つの FASTQ に混在）を
[xengsort](https://gitlab.com/genomeinformatics/xengsort) で host / graft に
分類し、pic のマッピング以降のパイプラインに接続するための手順です。

`xengsort` / `ntcard` / `gffread` は pic の conda パッケージに同梱されるため、
通常は追加インストール不要です。

## 全体像

```text
demux 済み FASTQ
      │  pic xenograft --index <name>
      ▼
xengsort classify (host / graft / both / neither / ambiguous)
      │  規約: graft -> graft genome, host -> host genome
      ▼
pic mapping --demux-fastq-dir <out>/classified/<fraction>  (画分ごと)
      ▼
pic deseq2 / enrich / report   (通常どおり)
```

## 1. xengsort index の作成（初回のみ）

host（-H）と graft（-G）のゲノム（必要なら cDNA も）から index を作成・登録する。
`-n`（k-mer 数）は省略時 `ntcard` の F0 から自動推定する。

```bash
pic build-xenograft-index \
  --name human_on_mouse \
  --host-genome mm10  --graft-genome hg38 \
  --host-fasta  mm10.fa.gz  --host-cdna  mm10.cdna.fa.gz \
  --graft-fasta hg38.fa.gz  --graft-cdna hg38.cdna.fa.gz
```

- 出力: `<PIC_LIB>/xengsort_index/<name>/<name>.hash|.info`
- 登録: `<PIC_LIB>/register/xengsort_index_map.tsv`
- 登録済み index は `pic build-xenograft-index --help` / `pic xenograft --help` の末尾に一覧表示される。

cDNA FASTA が無い生物種は、`--host-gtf` / `--graft-gtf` に GTF を渡すと
`gffread` でゲノム FASTA から cDNA を自動生成して index に含める:

```bash
pic build-xenograft-index \
  --name human_on_marmo \
  --host-genome calJac4 --graft-genome hg38 \
  --host-fasta  calJac4.fa.gz --host-gtf  calJac4.gtf.gz \
  --graft-fasta hg38.fa.gz    --graft-cdna hg38.cdna.fa.gz
```

## 2. 分類

```bash
pic xenograft \
  --index human_on_mouse \
  --demux-fastq-dir /path/demux \
  --out-dir /path/out \
  --fraction both        # graft / host / both (default: both)
```

出力:

- `<out>/xengsort/<sample>-{host,graft,both,neither,ambiguous}.fq.gz`
- `<out>/classified/<fraction>/<sample>.fastq.gz` … マッピング入力用
- `<out>/xenograft_classify_summary__<run>.tsv` … 画分ごとのリード数（QC）

## 3. マッピング以降（画分ごとに genome を指定）

規約どおり **graft は graft genome、host は host genome** にマッピングする。

```bash
# graft 画分 -> graft genome (例 hg38)
pic mapping --demux-fastq-dir /path/out/classified/graft \
  --sample-sheet <graft 用サンプルシート (genome=hg38)> --out-dir /path/out/graft_map

# host 画分 -> host genome (例 mm10)
pic mapping --demux-fastq-dir /path/out/classified/host \
  --sample-sheet <host 用サンプルシート (genome=mm10)> --out-dir /path/out/host_map
```

以降は通常どおり `pic deseq2` / `pic enrich` / `pic report`（または `pic all`）を実行する。

## 依存

pic conda パッケージ同梱: `xengsort` / `ntcard` / `gffread` / `pigz`。
