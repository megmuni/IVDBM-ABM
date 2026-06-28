#!/bin/bash
# Slurm batch script for ctest — submit via scripts/submit_tests.sh, not directly.
# Required env (set by submit_tests.sh): REPO_ROOT, BUILD_DIR, CTEST_ARGS
# Optional env: TEST_OUTPUT_DIR, IVDBM_PROFILE, TEST_SUITE

set -uo pipefail

fail() {
  echo "tests.sbatch: $*" >&2
  exit 1
}

cd "${REPO_ROOT:?Set REPO_ROOT when submitting this job.}" || fail "cannot cd to REPO_ROOT=${REPO_ROOT}"

TEST_DEFAULTS="${REPO_ROOT}/scripts/test_defaults.env"
[[ -f "${TEST_DEFAULTS}" ]] || fail "missing ${TEST_DEFAULTS}"
# shellcheck source=test_defaults.env
source "${TEST_DEFAULTS}"

BUILD_DIR="${BUILD_DIR:-${IVDBM_TEST_BUILD_DIR}}"
BUILD_DIR="${BUILD_DIR:-build}"
if [[ "${BUILD_DIR}" != /* ]]; then
  BUILD_DIR="${REPO_ROOT}/${BUILD_DIR}"
fi

TEST_OUTPUT_DIR="${TEST_OUTPUT_DIR:-output/tests_${SLURM_JOB_ID}}"
if [[ "${TEST_OUTPUT_DIR}" != /* ]]; then
  TEST_OUTPUT_DIR="${REPO_ROOT}/${TEST_OUTPUT_DIR}"
fi
CTEST_LOG="${TEST_OUTPUT_DIR}/ctest.log"
CTEST_ARGS="${CTEST_ARGS:-}"
TEST_SUITE="${TEST_SUITE:-${IVDBM_TEST_DEFAULT_SUITE}}"

mkdir -p "${TEST_OUTPUT_DIR}" logs || fail "cannot create ${TEST_OUTPUT_DIR} or logs/"

[[ -d "${BUILD_DIR}" ]] || fail "build directory missing: ${BUILD_DIR}"
[[ -f "${BUILD_DIR}/CTestTestfile.cmake" ]] || \
  fail "no CTest config in ${BUILD_DIR} — run ./scripts/build_drac.sh --tests"

TEMPLATE="${REPO_ROOT}/configFiles/chemical_environment.template.json"
[[ -f "${TEMPLATE}" ]] || \
  fail "missing ${TEMPLATE} (required by chemistry/world tests)"

export LD_LIBRARY_PATH="${BUILD_DIR}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

if [[ -f /etc/profile.d/modules.sh ]]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
fi
module load StdEnv/2023 2>&1 || true
if [[ "${IVDBM_PROFILE:-cpu}" == "gpu" ]]; then
  module load cuda 2>&1 || module load cuda/12.2 2>&1 || true
fi

echo "Test job ${SLURM_JOB_ID:-?} starting on $(hostname)"
echo "  REPO_ROOT=${REPO_ROOT}"
echo "  BUILD_DIR=${BUILD_DIR}"
echo "  TEST_OUTPUT_DIR=${TEST_OUTPUT_DIR}"
echo "  TEST_SUITE=${TEST_SUITE}"
echo "  IVDBM_PROFILE=${IVDBM_PROFILE:-cpu}"
echo "  CTEST_ARGS=${CTEST_ARGS:-<none>}"

set +e
export CTEST_OUTPUT_ON_FAILURE=1
# shellcheck disable=SC2086
ctest --test-dir "${BUILD_DIR}" ${CTEST_ARGS} 2>&1 | tee "${CTEST_LOG}"
status=${PIPESTATUS[0]}
set -e

echo "ctest finished with exit code ${status}"
exit "${status}"
