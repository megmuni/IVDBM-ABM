#!/usr/bin/env bash
# Configure and build testRun on a DRAC (Alliance Canada) login node.
#
# Usage:
#   ./scripts/build_drac.sh              # CPU-only (no cuda module required)
#   ./scripts/build_drac.sh --cuda       # GPU diffusion (loads cuda module)
#   ./scripts/build_drac.sh --cuda --arch 90

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build"
BUILD_TYPE="Release"
ENABLE_CUDA=0
CUDA_ARCH="90"
BUILD_JOBS="${BUILD_JOBS:-$(nproc 2>/dev/null || echo 4)}"

usage() {
  cat <<'EOF'
Build IVDBM-ABM on a DRAC login node.

Usage:
  build_drac.sh [--cuda] [--arch SM] [--debug] [--clean]

Options:
  --cuda        Enable GPU diffusion (-DDIFFUSION3D_CUDA=ON; loads cuda module)
  --arch SM     CUDA architecture for --cuda (default: 90 for H100)
  --debug       Use Debug build type instead of Release
  --clean       Remove build/ before configuring
  -h, --help    Show this help

CPU-only example (login node, no nvcc needed):
  ./scripts/build_drac.sh

GPU example (must load cuda before cmake):
  ./scripts/build_drac.sh --cuda --arch 90
EOF
}

die() {
  echo "build_drac.sh: $*" >&2
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --cuda)
      ENABLE_CUDA=1
      shift
      ;;
    --arch)
      [[ $# -ge 2 ]] || die "missing value for $1"
      CUDA_ARCH="$2"
      shift 2
      ;;
    --debug)
      BUILD_TYPE="Debug"
      shift
      ;;
    --clean)
      rm -rf "${BUILD_DIR}"
      shift
      ;;
    *)
      die "unknown option: $1 (try --help)"
      ;;
  esac
done

module load StdEnv/2023 2>/dev/null || true
module load gcc/12.3 2>/dev/null || module load gcc 2>/dev/null || true
module load cmake 2>/dev/null || true

if [[ "${ENABLE_CUDA}" -eq 1 ]]; then
  module load cuda/12.2 2>/dev/null || module load cuda 2>/dev/null || \
    die "could not load a cuda module; run 'module avail cuda' and load one manually"
  command -v nvcc >/dev/null 2>&1 || die "nvcc not on PATH after loading cuda module"
fi

mkdir -p "${BUILD_DIR}"
cd "${REPO_ROOT}"

CMAKE_ARGS=(
  -S .
  -B "${BUILD_DIR}"
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
  -DBUILD_SRC_TESTS=OFF
)

if [[ "${ENABLE_CUDA}" -eq 1 ]]; then
  CMAKE_ARGS+=(
    -DDIFFUSION3D_CUDA=ON
    -DDIFFUSION3D_CUDA_ARCHITECTURES="${CUDA_ARCH}"
  )
fi

echo "==> cmake ${CMAKE_ARGS[*]}"
cmake "${CMAKE_ARGS[@]}"

echo "==> cmake --build ${BUILD_DIR} --target testRun -j${BUILD_JOBS}"
cmake --build "${BUILD_DIR}" --target testRun -j"${BUILD_JOBS}"

echo "Built: ${BUILD_DIR}/bin/testRun"
