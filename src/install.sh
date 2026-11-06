#!/usr/bin/env bash
# pic の Rust のプログラム (pic-demux・pic-aggregate) をまとめて build して、<bin_dir> に置く。
#   bash src/install.sh [<bin_dir>]     (既定: この toolkit の bin/、つまり bin/pic の隣)
# bin/pic は自分の隣の bin/pic-demux・bin/pic-aggregate を最初に使う (無ければ PATH のもの)。
# toolkit の写し (git worktree など) ごとにこれを回せば、版ごとに持てる。conda パッケージ
# (src/recipe/build.sh) は <bin_dir> に "$PREFIX/bin" を渡す。
# どちらも中身の hash ごとに build したものを取っておくので、変わっていなければ写すだけ。
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bin_dir="${1:-$(cd "${here}/.." && pwd)/bin}"

bash "${here}/pic-demux/install.sh" "$bin_dir"
bash "${here}/pic-aggregate/install.sh" "$bin_dir"
