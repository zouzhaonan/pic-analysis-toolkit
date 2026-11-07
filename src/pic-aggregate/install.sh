#!/usr/bin/env bash
# pic-aggregate (pic aggregate のメタジーン集計のプログラム) を build して、<bin_dir> に置く。
# ふつうは src/install.sh (pic-demux と pic-aggregate をまとめて置く) から呼ぶ。
#   bash src/pic-aggregate/install.sh [<bin_dir>]     (既定: この toolkit の bin/、つまり bin/pic の隣)
# bin/pic は、自分の隣の bin/pic-aggregate を最初に使う (無ければ PATH の pic-aggregate)。toolkit の写し
# (git worktree など) ごとに install.sh を回せば、版ごとに pic-aggregate を持てる。bin/pic-aggregate は git に入れない。
# conda パッケージ (src/recipe/build.sh) は <bin_dir> に "$PREFIX/bin" を渡す。
# 要るもの: cargo (rust) と、リンクに C コンパイラ (cc か $CC)。C のライブラリは使わない (zlib は zlib-rs)。
# 置かなくても pic aggregate は動く (無ければ bigWigAverageOverBed の道に戻る。結果は同じ)。
#
# 速くするため、build したものを src/pic-aggregate の中身 (Cargo.toml・Cargo.lock・src/) の hash ごとに
# ${PIC_AGGREGATE_CACHE:-${XDG_CACHE_HOME:-~/.cache}/pic-aggregate}/<hash>/pic-aggregate に取っておく。同じ中身なら
# build せずに写すだけ (cargo も要らない)。build の作業場所 (<cache>/target) も写しどうしで共有する。
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bin_dir="${1:-$(cd "${here}/../.." && pwd)/bin}"
cache_root="${PIC_AGGREGATE_CACHE:-${XDG_CACHE_HOME:-${HOME}/.cache}/pic-aggregate}"

sha256() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum; else shasum -a 256; fi
}
src_hash="$(
  cd "$here"
  { printf '%s\n' Cargo.toml Cargo.lock; find src -type f | LC_ALL=C sort; } |
    while IFS= read -r f; do printf '%s\0' "$f"; cat "$f"; done | sha256
)"
src_hash="${src_hash%% *}"
src_hash="${src_hash:0:16}"
cached="${cache_root}/${src_hash}/pic-aggregate"

if [[ ! -x "$cached" ]]; then
  if ! command -v cargo >/dev/null 2>&1; then
    echo "[ERROR] cargo (rust) が見つかりません。環境 pic に入れてください: micromamba install -n pic -c conda-forge rust" >&2
    exit 1
  fi
  # cc が無く、conda の C コンパイラ ($CC、環境 pic を有効にすると決まる) だけがあるとき (workstation)、
  # それを Rust のリンカーにも使う。
  if ! command -v cc >/dev/null 2>&1 && [[ -n "${CC:-}" ]]; then
    export RUSTFLAGS="${RUSTFLAGS:-} -C linker=${CC}"
  fi
  target_dir="${CARGO_TARGET_DIR:-${cache_root}/target}"
  mkdir -p "${cache_root}/${src_hash}"
  # 同じ target を同時に使うときは cargo が自分で待つ。
  cargo build --release --locked --manifest-path "${here}/Cargo.toml" --target-dir "$target_dir"
  # conda-forge の rust の activation は CARGO_BUILD_TARGET を決めることがあり、そのときは target/<triple>/release に出る。
  built="${target_dir}/${CARGO_BUILD_TARGET:+${CARGO_BUILD_TARGET}/}release/pic-aggregate"
  install -m 755 "$built" "${cached}.tmp.$$"
  mv "${cached}.tmp.$$" "$cached"
fi

mkdir -p "$bin_dir"
install -m 755 "$cached" "${bin_dir}/pic-aggregate.tmp.$$"
mv "${bin_dir}/pic-aggregate.tmp.$$" "${bin_dir}/pic-aggregate"
"${bin_dir}/pic-aggregate" --version
