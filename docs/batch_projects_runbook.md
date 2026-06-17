# プロジェクト一括処理ランブック (Claude 向け手順書)

このドキュメントは、`Book2.txt` のような **6 列テーブル** を渡されたときに Claude が
プロジェクト単位で `pic all` を実行するための手順書です。`pic all` 自体は単一サンプル
シート専用のシンプルなコマンドなので、複数プロジェクトへの分割・genome 整備の判断・
ユーザーへの許可取りはこの手順に従って Claude が行います。

## 入力

### 6 列テーブル (ヘッダ無し・タブ区切り)

| col | 意味 |
|-----|------|
| 1 | fastq_prefix |
| 2 | genome |
| 3 | barcode |
| 4 | sample |
| 5 | group |
| 6 | project (出力フォルダ名) |

1〜5 列目は `pic mapping` のサンプルシートと同形式。例 (`/Users/zou/Downloads/Book2.txt`):

```
1_S1   xenbase_v10.1   CATGAG   GFP_posi_1   GFP_posi   pic163
3_S3   mm10            TCGAAG   WT_Contact   WT_Contact pic164
...
```

### ユーザーから別途受け取るもの

- `raw-fastq-dir`: 生ペア FASTQ (R1/R2) の置き場所。
- `out-dir-base`: 出力のベースディレクトリ (この下に project ごとのフォルダを作る)。

これらが未指定なら **実行前に必ずユーザーへ確認する**。

## 手順

### 1. プロジェクト毎に分割

col6 (project) でグループ化し、各プロジェクト `P` について 1〜5 列目を取り出して
`<out-dir-base>/<P>/sample_sheet.tsv` を書き出す。先頭に必ずヘッダ行を付ける:

```
fastq_prefix	genome	barcode	sample	group
```

### 2. genome 整備チェック

各プロジェクトに含まれる genome について、pic の参照 lib が揃っているか確認する。
lib のルートは `${PIC_LIB:-$HOME/local/lib/pic}`。

| ステップ | 必要なもの | 確認先 |
|----------|-----------|--------|
| mapping | genome_map 登録 + index/gtf/chrom | `register/genome_map.tsv` に行があり、`hisat2_index/<g>/genome.1.ht2` (または `.1.ht2l`) / `gtf/<g>.gtf`(`.gz` 可) / `chrom_size/<g>.chrom.sizes` が存在 |
| deseq2 | biomart 登録 | `register/biomart_lookup.tsv` に `<g>` の行 |
| enrich | enrichment マップ + lib | `register/genome_enrichment_map.tsv` に `<g>` + enrichment lib 準備済み |

簡易確認コマンド例:

```bash
LIB="${PIC_LIB:-$HOME/local/lib/pic}"
pic mapping --help          # 末尾に登録済み pic genome 一覧が出る
pic manage-biomart --help   # 末尾に登録済み biomart genome 一覧が出る
ls "$LIB/hisat2_index/<g>/" "$LIB/register/"
```

- **全て揃っている → そのまま続行。**
- **欠けている → 整備コマンドを提示し、`AskUserQuestion` でユーザーの許可を得てから実行する。** 許可後に再チェックして揃ったことを確認する。

  | 欠けているもの | 整備コマンド |
  |----------------|--------------|
  | mapping 参照 | `pic build-genome --genome <g> --fasta <path> --gtf <path>` |
  | biomart | `pic manage-biomart --register --genome <g> --dataset <dataset> [--source <src>]`(まず `pic manage-biomart --list-datasets` で dataset 名を確認) |
  | enrichment lib | `pic enrich --prepare-libs --genome <g>` |

  FASTA/GTF や dataset 名など Claude が知り得ない入力はユーザーに尋ねる。

### 3. 非標準 genome 名は必ずユーザーに確認

`Human/Mouse` のようにスラッシュ・空白・複数種混成を思わせる名前は、そのまま lib 名
として使えない(混成参照か別名の可能性)。**該当時は勝手に置換・推測せず、登録名や参照
ファイルの方針をユーザーに確認してから進める。**(ユーザー方針: 非典型的な genome 名が
あったら都度確認する。)

### 4. 実行

各プロジェクトについて:

```bash
pic all \
  --sample-sheet <out-dir-base>/<P>/sample_sheet.tsv \
  --raw-fastq-dir <raw-fastq-dir> \
  --run-name <P> \
  --out-dir <out-dir-base>/<P>
```

出力:

- `<out-dir-base>/<P>/` … mapping 成果物 (deftable, counts, bam, bw …)
- `<out-dir-base>/<P>/deseq2/<genome>/` … DESeq2 結果
- `<out-dir-base>/<P>/enrich/<genome>/` … ORA/GSEA 結果

`pic all` は `set -euo pipefail` で動くため、どこかで失敗すればそのプロジェクトは
そこで停止する。失敗したら原因(多くは genome 未整備)を確認して整備し直す。

### 5. サマリ報告

完了後、プロジェクト × genome の成否、整備が必要だった genome、要確認のためスキップ
した項目を表でユーザーに報告する。

## 参考: Book2.txt のプロジェクト/genome 対応(2026 時点の例)

| project | genome |
|---------|--------|
| pic163 | xenbase_v10.1 |
| pic164 | mm10 |
| DrSawada | Human/Mouse(※非標準名 → 要確認) |
| pic165 | mat1.0 |
| pic169 | mm10 |
| pic170 | mm10 |
| pic171 | mm10 |
| dev026 | hg38 |
