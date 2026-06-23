#!/usr/bin/env bash
# Submit IVDBM-ABM Catch2/CTest jobs to DRAC (Slurm).
#
# Usage:
#   ./scripts/submit_tests.sh $EMAIL --suite chemistry
#   ./scripts/submit_tests.sh $EMAIL --suite gpu --profile gpu

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SBATCH_SCRIPT="${SCRIPT_DIR}/tests.sbatch"
TEST_DEFAULTS="${SCRIPT_DIR}/test_defaults.env"

[[ -f "${TEST_DEFAULTS}" ]] || {
  echo "submit_tests.sh: missing ${TEST_DEFAULTS}" >&2
  exit 1
}
# shellcheck source=test_defaults.env
source "${TEST_DEFAULTS}"

EMAIL=""
ACCOUNT="${IVDBM_TEST_DEFAULT_ACCOUNT}"
ARRAY="${IVDBM_TEST_DEFAULT_ARRAY}"
TIME="${IVDBM_TEST_DEFAULT_TIME}"
CPUS="${IVDBM_TEST_DEFAULT_CPUS}"
GPUS=""
MEM="${IVDBM_TEST_DEFAULT_MEM}"
PROFILE="${IVDBM_TEST_DEFAULT_PROFILE}"
SUITE="${IVDBM_TEST_DEFAULT_SUITE}"
BUILD_DIR="${IVDBM_TEST_BUILD_DIR}"
OUTPUT_DIR=""
CTEST_EXTRA=()
DRY_RUN=0
GPUS_EXPLICIT=0

usage() {
  cat <<'EOF'
Submit Catch2/CTest jobs to DRAC (Slurm).

Usage:
  submit_tests.sh EMAIL [options]
  submit_tests.sh --email EMAIL [options]

Required:
  EMAIL                 Address for Slurm mail notifications (--mail-user)

Options:
  --suite NAME          Test suite: chemistry, diffusion3d, world, cpu, gpu, all
                        (default: cpu)
  --profile cpu|gpu     Slurm profile (default: cpu; use gpu for --suite gpu)
  --build-dir PATH      CMake build directory (default: build)
  --output-dir PATH     Test artifacts (default: output/tests_<SLURM_JOB_ID>)
  --account ACCOUNT     Slurm allocation (default: def-nicoleli)
  --array RANGE         Job array range (default: 0-0)
  --time D-HH:MM:SS     Wall time (default: 0-00:30:00)
  --cpus N              CPUs per task (default: 8)
  --gpus SPEC           GPUs per node (overrides --profile)
  --mem SIZE            Memory per node (default: 16000M)
  --dry-run             Print sbatch command without submitting
  -h, --help            Show this help

Build tests first:
  ./scripts/build_drac.sh --tests              # Catch2 v2 module or FetchContent fallback
  ./scripts/build_drac.sh --cuda --tests       # include GPU tests

Examples:
  ./scripts/submit_tests.sh $EMAIL --suite chemistry
  ./scripts/submit_tests.sh $EMAIL --suite cpu
  ./scripts/submit_tests.sh $EMAIL --suite gpu --profile gpu
  ./scripts/submit_tests.sh $EMAIL --suite all --time 1-00:00:00
  ./scripts/submit_tests.sh $EMAIL -- --output-on-failure -V
EOF
}

die() {
  echo "submit_tests.sh: $*" >&2
  exit 1
}

warn() {
  echo "submit_tests.sh: warning: $*" >&2
}

valid_email() {
  [[ "$1" =~ ^[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}$ ]]
}

apply_profile() {
  case "${PROFILE}" in
    cpu)
      GPUS="${IVDBM_PROFILE_CPU_GPUS}"
      ;;
    gpu)
      GPUS="${IVDBM_PROFILE_GPU_GPUS}"
      ;;
    *)
      die "unknown profile: ${PROFILE} (use cpu or gpu)"
      ;;
  esac
}

resolve_build_dir() {
  local path="$1"
  if [[ "${path}" != /* ]]; then
    path="${REPO_ROOT}/${path}"
  fi
  printf '%s' "${path}"
}

build_has_cuda() {
  [[ -f "${BUILD_DIR}/CMakeCache.txt" ]] \
    && grep -q '^DIFFUSION3D_CUDA:BOOL=ON' "${BUILD_DIR}/CMakeCache.txt"
}

build_has_tests() {
  [[ -f "${BUILD_DIR}/CTestTestfile.cmake" ]]
}

suite_to_ctest_args() {
  case "${SUITE}" in
    chemistry)
      CTEST_ARGS="-L chemistry"
      ;;
    diffusion3d)
      CTEST_ARGS="-L diffusion3d"
      ;;
    world)
      CTEST_ARGS="-L world"
      ;;
    cpu)
      CTEST_ARGS="-L cpu"
      ;;
    gpu)
      CTEST_ARGS="-L gpu"
      ;;
    all)
      CTEST_ARGS=""
      ;;
    *)
      die "unknown suite: ${SUITE} (use chemistry, diffusion3d, world, cpu, gpu, or all)"
      ;;
  esac
}

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    -e|--email)
      [[ $# -ge 2 ]] || die "missing value for $1"
      EMAIL="$2"
      shift 2
      ;;
    --suite)
      [[ $# -ge 2 ]] || die "missing value for $1"
      SUITE="$2"
      shift 2
      ;;
    --profile)
      [[ $# -ge 2 ]] || die "missing value for $1"
      PROFILE="$2"
      shift 2
      ;;
    --build-dir)
      [[ $# -ge 2 ]] || die "missing value for $1"
      BUILD_DIR="$2"
      shift 2
      ;;
    --output-dir)
      [[ $# -ge 2 ]] || die "missing value for $1"
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --account)
      [[ $# -ge 2 ]] || die "missing value for $1"
      ACCOUNT="$2"
      shift 2
      ;;
    --array)
      [[ $# -ge 2 ]] || die "missing value for $1"
      ARRAY="$2"
      shift 2
      ;;
    --time)
      [[ $# -ge 2 ]] || die "missing value for $1"
      TIME="$2"
      shift 2
      ;;
    --cpus)
      [[ $# -ge 2 ]] || die "missing value for $1"
      CPUS="$2"
      shift 2
      ;;
    --gpus)
      [[ $# -ge 2 ]] || die "missing value for $1"
      GPUS="$2"
      GPUS_EXPLICIT=1
      shift 2
      ;;
    --mem)
      [[ $# -ge 2 ]] || die "missing value for $1"
      MEM="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --)
      shift
      while [[ $# -gt 0 ]]; do
        CTEST_EXTRA+=("$1")
        shift
      done
      break
      ;;
    -*)
      die "unknown option: $1 (try --help)"
      ;;
    *)
      if [[ -n "${EMAIL}" ]]; then
        die "unexpected argument: $1 (email already set to ${EMAIL})"
      fi
      EMAIL="$1"
      shift
      ;;
  esac
done

[[ -n "${EMAIL}" ]] || die "EMAIL is required (positional or --email)"
valid_email "${EMAIL}" || die "invalid email address: ${EMAIL}"

BUILD_DIR="$(resolve_build_dir "${BUILD_DIR}")"

if [[ "${SUITE}" == "gpu" && "${PROFILE}" == "cpu" ]]; then
  PROFILE="gpu"
fi

if [[ "${GPUS_EXPLICIT}" -eq 0 ]]; then
  apply_profile
fi

suite_to_ctest_args

if [[ ${#CTEST_EXTRA[@]} -gt 0 ]]; then
  CTEST_ARGS="${CTEST_ARGS} ${CTEST_EXTRA[*]}"
fi

[[ -f "${SBATCH_SCRIPT}" ]] || die "missing batch script: ${SBATCH_SCRIPT}"

if ! build_has_tests; then
  die "tests not configured in ${BUILD_DIR} — run: ./scripts/build_drac.sh --tests"
fi

if [[ "${SUITE}" == "gpu" || "${SUITE}" == "all" ]]; then
  if ! build_has_cuda; then
    warn "suite=${SUITE} may include GPU tests but build was not configured with -DDIFFUSION3D_CUDA=ON"
    warn "rebuild with: ./scripts/build_drac.sh --cuda --tests"
  fi
fi

if [[ "${SUITE}" == "gpu" && "${PROFILE}" != "gpu" ]]; then
  warn "suite=gpu usually needs --profile gpu to request a GPU node"
fi

mkdir -p "${REPO_ROOT}/logs"

EXPORT_VARS="NONE,REPO_ROOT=${REPO_ROOT},BUILD_DIR=${BUILD_DIR},CTEST_ARGS=${CTEST_ARGS},IVDBM_PROFILE=${PROFILE},TEST_SUITE=${SUITE}"
if [[ -n "${OUTPUT_DIR}" ]]; then
  EXPORT_VARS="${EXPORT_VARS},TEST_OUTPUT_DIR=${OUTPUT_DIR}"
fi

SBATCH_ARGS=(
  --job-name=ivdbm-tests
  --account="${ACCOUNT}"
  --mail-user="${EMAIL}"
  --mail-type=ALL
  --time="${TIME}"
  --nodes=1
  --cpus-per-task="${CPUS}"
  --mem="${MEM}"
  --array="${ARRAY}"
  --chdir="${REPO_ROOT}"
  --output="${REPO_ROOT}/logs/tests_%a_%j.out"
  --error="${REPO_ROOT}/logs/tests_%a_%j.err"
  --export="${EXPORT_VARS}"
)

if [[ -n "${GPUS}" ]]; then
  SBATCH_ARGS+=(--gpus-per-node="${GPUS}")
fi

SBATCH_ARGS+=("${SBATCH_SCRIPT}")

if [[ "${DRY_RUN}" -eq 1 ]]; then
  printf 'Would run: sbatch'
  printf ' %q' "${SBATCH_ARGS[@]}"
  printf '\n'
  echo "  suite=${SUITE}  ctest args: ${CTEST_ARGS:-<all>}"
  exit 0
fi

command -v sbatch >/dev/null 2>&1 || die "sbatch not found (log in to a DRAC cluster login node first)"

SUBMIT_OUTPUT="$(sbatch "${SBATCH_ARGS[@]}")"
echo "${SUBMIT_OUTPUT}"
JOB_ID="${SUBMIT_OUTPUT##* }"
echo "Build dir: ${BUILD_DIR}"
echo "Suite: ${SUITE} (${CTEST_ARGS:-all tests})"
echo "Profile: ${PROFILE}${GPUS:+ (${GPUS})}"
if [[ -n "${OUTPUT_DIR}" ]]; then
  echo "Output dir (override): ${OUTPUT_DIR}"
else
  echo "Output dir: output/tests_<SLURM_JOB_ID>"
fi
echo "Logs: logs/tests_<array>_${JOB_ID}.out"
echo "Mail notifications for ${EMAIL} will be sent by Slurm for job ${JOB_ID}."
