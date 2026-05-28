#!/usr/bin/env bash
set -euo pipefail

pkg_root="${SRC_DIR}"
install_root="${PREFIX}/pic"

mkdir -p "${PREFIX}/bin"
mkdir -p "${install_root}"

cp -R "${pkg_root}/bin" "${install_root}/"
cp -R "${pkg_root}/scripts" "${install_root}/"
cp -R "${pkg_root}/help" "${install_root}/"

cat > "${PREFIX}/bin/pic" <<'WRAP'
#!/usr/bin/env bash
set -euo pipefail

PIC_ROOT="${CONDA_PREFIX}/pic"
exec "${PIC_ROOT}/bin/pic" "$@"
WRAP

chmod +x "${PREFIX}/bin/pic"
chmod +x "${install_root}/bin/pic"
