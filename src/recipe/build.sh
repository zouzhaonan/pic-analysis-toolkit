#!/usr/bin/env bash
set -euo pipefail

pkg_root="${SRC_DIR}"
install_root="${PREFIX}/pic"

mkdir -p "${PREFIX}/bin"
mkdir -p "${install_root}"

cp -R "${pkg_root}/bin" "${install_root}/"
cp -R "${pkg_root}/scripts" "${install_root}/"
cp -R "${pkg_root}/help" "${install_root}/"
cp "${pkg_root}/VERSION" "${install_root}/"

# pic demux の振り分けのプログラム (Rust)。環境の bin に置く。
CARGO_TARGET_DIR="${SRC_DIR}/target-pic-demux" PIC_DEMUX_CACHE="${SRC_DIR}/cache-pic-demux" bash "${pkg_root}/src/pic-demux/install.sh" "${PREFIX}/bin"

cat > "${PREFIX}/bin/pic" <<'WRAP'
#!/usr/bin/env bash
set -euo pipefail

PIC_ROOT="${CONDA_PREFIX}/pic"
exec "${PIC_ROOT}/bin/pic" "$@"
WRAP

chmod +x "${PREFIX}/bin/pic"
chmod +x "${install_root}/bin/pic"
