#!/usr/bin/env bash
# pic-aggregate (pic aggregate のメタジーン集計のプログラム) を build して、<bin_dir> に置く。
#   bash src/pic-aggregate/install.sh [<bin_dir>]     (既定: "$CONDA_PREFIX/bin"、つまり環境 pic の bin)
# 要るもの: cargo (rust)。C のライブラリは使わない (zlib は Rust の zlib-rs)。リンクに C コンパイラ (cc か $CC) が要る。
# build の作業場所は src/pic-aggregate/target (git に入れない)。変わっていなければ build はすぐ終わる。
# 置かなくても pic aggregate は動く (PATH に無ければ bigWigAverageOverBed の道に戻る。結果は同じ)。
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bin_dir="${1:-${CONDA_PREFIX:?CONDA_PREFIX is not set; activate the pic environment or give <bin_dir>}/bin}"
target_dir="${CARGO_TARGET_DIR:-${here}/target}"

# cc が無く、conda の C コンパイラ ($CC、環境 pic を有効にすると決まる) だけがあるとき (workstation)、
# それを Rust のリンカーにも使う。
if ! command -v cc >/dev/null 2>&1 && [[ -n "${CC:-}" ]]; then
  export RUSTFLAGS="${RUSTFLAGS:-} -C linker=${CC}"
fi

cargo build --release --locked --manifest-path "${here}/Cargo.toml" --target-dir "$target_dir"

# conda-forge の rust の activation は CARGO_BUILD_TARGET を決めることがあり、そのときは target/<triple>/release に出る。
built="${target_dir}/${CARGO_BUILD_TARGET:+${CARGO_BUILD_TARGET}/}release/pic-aggregate"
mkdir -p "$bin_dir"
install -m 755 "$built" "${bin_dir}/pic-aggregate.tmp"
mv "${bin_dir}/pic-aggregate.tmp" "${bin_dir}/pic-aggregate"
"${bin_dir}/pic-aggregate" --version
