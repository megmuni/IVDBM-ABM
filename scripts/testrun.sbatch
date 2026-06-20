#!/bin/bash
# Slurm batch script for testRun — submit via scripts/submit_testrun.sh, not directly.
# Required env (set by submit_testrun.sh): REPO_ROOT, TESTRUN
# Optional env: OUTPUT_DIR, NUMTICKS, WXW, WYW, WZW, INPUTFILE, OUTPUTFILE, CHEM_CONFIG

set -uo pipefail

fail() {
  echo "testrun.sbatch: $*" >&2
  exit 1
}

cd "${REPO_ROOT:?Set REPO_ROOT when submitting this job.}" || fail "cannot cd to REPO_ROOT=${REPO_ROOT}"

OUTPUT_DIR="${OUTPUT_DIR:-output/job_${SLURM_JOB_ID}}"
INPUTFILE="${INPUTFILE:-configFiles/config_scaffold.txt}"
CHEM_CONFIG="${CHEM_CONFIG:-configFiles/chemical_environment.json}"
NUMTICKS="${NUMTICKS:-24}"
WXW="${WXW:-1}"
WYW="${WYW:-1}"
WZW="${WZW:-1}"
OUTPUTFILE="${OUTPUTFILE:-${OUTPUT_DIR}/Output_Biomarkers.csv}"

mkdir -p "${OUTPUT_DIR}" logs || fail "cannot create ${OUTPUT_DIR} or logs/"

[[ -x "${TESTRUN:?Set TESTRUN when submitting this job.}" ]] \
  || fail "testRun not executable: ${TESTRUN}"
[[ -f "${INPUTFILE}" ]] || fail "input file missing: ${REPO_ROOT}/${INPUTFILE}"
[[ -f "${CHEM_CONFIG}" ]] \
  || fail "chemical config missing: ${REPO_ROOT}/${CHEM_CONFIG} (copy configFiles/chemical_environment.template.json)"

export LD_LIBRARY_PATH="${REPO_ROOT}/build/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
export IVDBM_OUTPUT_DIR="${OUTPUT_DIR}"
export IVDBM_RUN_PARAMS_JSON="${OUTPUT_DIR}/run_params.json"

if [[ -f /etc/profile.d/modules.sh ]]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
fi
module load StdEnv/2023 2>&1 || true
module load cuda 2>&1 || module load cuda/12.2 2>&1 || true

echo "Job ${SLURM_JOB_ID:-?} starting on $(hostname)"
echo "  REPO_ROOT=${REPO_ROOT}"
echo "  TESTRUN=${TESTRUN}"
echo "  OUTPUT_DIR=${OUTPUT_DIR}"
echo "  INPUTFILE=${INPUTFILE}"
echo "  CHEM_CONFIG=${CHEM_CONFIG}"
echo "  NUMTICKS=${NUMTICKS}  OMP_NUM_THREADS=${OMP_NUM_THREADS}"

set -e
"${TESTRUN}" \
  --output-dir "${OUTPUT_DIR}" \
  --numticks "${NUMTICKS}" \
  --inputfile "${INPUTFILE}" \
  --chem-config "${CHEM_CONFIG}" \
  --outputfile "${OUTPUTFILE}" \
  --wxw "${WXW}" \
  --wyw "${WYW}" \
  --wzw "${WZW}"
status=$?
echo "testRun finished with exit code ${status}"
exit "${status}"
