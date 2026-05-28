# PIC Post-Install Setup

`pic` package includes `lib/help` only.
Runtime data (genome index/cache/register/enrichment) must be prepared externally and mounted.

## 1. Prepare runtime lib directory
Prepare a runtime lib directory path.  
An empty directory is acceptable; required subdirectories are created as needed during execution.

## 2. Runtime lib path
- Default path: `$HOME/local/lib/pic` (auto-set by `pic` when `PIC_LIB` is unset)
- Override only when needed:
```bash
export PIC_LIB=/path/to/mounted/lib
```

---

# PIC インストール後設定

`pic` パッケージには `lib/help` のみ含まれます。  
実行時データ（genome index/cache/register/enrichment）は外部で準備してマウントする必要があります。

## 1. 実行時 lib ディレクトリを準備
実行時に使う `lib` ディレクトリのパスを用意してください。  
空ディレクトリでも問題ありません。必要なサブディレクトリは実行時に作成されます。

## 2. 実行時 lib パス
- デフォルト: `$HOME/local/lib/pic`（`PIC_LIB` 未設定時は `pic` が自動設定）
- 必要な場合のみ上書き:
```bash
export PIC_LIB=/path/to/mounted/lib
```

## 3. 設定確認
`help` は lib なしでも表示できますが、解析コマンドには lib が必要です。

```bash
pic help
pic mapping --help
pic deseq2 --help
```

## 4. 永続化（任意）
シェル設定ファイル（`~/.bashrc` または `~/.zshrc`）に追加:

```bash
export PIC_LIB=/path/to/mounted/lib
```

## 3. Verify configuration
`help` works without runtime lib, but analysis commands require it.

```bash
pic help
pic mapping --help
pic deseq2 --help
```

## 4. Persistent setting (optional)
Add to your shell profile (`~/.bashrc` or `~/.zshrc`):

```bash
export PIC_LIB=/path/to/mounted/lib
```
