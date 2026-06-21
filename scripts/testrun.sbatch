#!/bin/bash
# Slurm batch script for testRun — submit via scripts/submit_testrun.sh, not directly.
# Required env (set by submit_testrun.sh): REPO_ROOT, TESTRUN
# Optional env: OUTPUT_DIR, NUMTICKS, WXW, WYW, WZW, INPUTFILE, OUTPUTFILE, CHEM_CONFIG, IVDBM_PROFILE

set -uo pipefail

fail() {
  echo "testrun.sbatch: $*" >&2
  exit 1
}

cd "${REPO_ROOT:?Set REPO_ROOT when submitting this job.}" || fail "cannot cd to REPO_ROOT=${REPO_ROOT}"

JOB_DEFAULTS="${REPO_ROOT}/scripts/job_defaults.env"
[[ -f "${JOB_DEFAULTS}" ]] || fail "missing ${JOB_DEFAULTS}"
# shellcheck source=job_defaults.env
source "${JOB_DEFAULTS}"

OUTPUT_DIR="${OUTPUT_DIR:-output/job_${SLURM_JOB_ID}}"
OUTPUTFILE="${OUTPUTFILE:-${OUTPUT_DIR}/Output_Biomarkers.csv}"
FINALIZE_SCRIPT="${REPO_ROOT}/scripts/finalize_run_params.py"
RUN_TIMING_JSON="${OUTPUT_DIR}/.run_timing.json"
TIME_FILE="${OUTPUT_DIR}/.testrun.time"

mkdir -p "${OUTPUT_DIR}" logs || fail "cannot create ${OUTPUT_DIR} or logs/"

[[ -x "${TESTRUN:?Set TESTRUN when submitting this job.}" ]] \
  || fail "testRun not executable: ${TESTRUN}"
[[ -f "${INPUTFILE}" ]] || fail "input file missing: ${REPO_ROOT}/${INPUTFILE}"
[[ -f "${CHEM_CONFIG}" ]] \
  || fail "chemical config missing: ${REPO_ROOT}/${CHEM_CONFIG} (copy configFiles/chemical_environment.template.json)"
[[ -f "${FINALIZE_SCRIPT}" ]] || fail "missing ${FINALIZE_SCRIPT}"

export LD_LIBRARY_PATH="${REPO_ROOT}/build/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
export IVDBM_OUTPUT_DIR="${OUTPUT_DIR}"
export IVDBM_RUN_PARAMS_JSON="${OUTPUT_DIR}/run_params.json"
export IVDBM_RUN_TIMING_JSON="${RUN_TIMING_JSON}"

if [[ -f /etc/profile.d/modules.sh ]]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
fi
module load StdEnv/2023 2>&1 || true
if [[ "${IVDBM_PROFILE:-cpu}" == "gpu" ]]; then
  module load cuda 2>&1 || module load cuda/12.2 2>&1 || true
fi

JOB_STARTED_AT="$(date -Iseconds)"
JOB_HOST="$(hostname)"

echo "Job ${SLURM_JOB_ID:-?} starting on ${JOB_HOST} at ${JOB_STARTED_AT}"
echo "  REPO_ROOT=${REPO_ROOT}"
echo "  TESTRUN=${TESTRUN}"
echo "  OUTPUT_DIR=${OUTPUT_DIR}"
echo "  INPUTFILE=${INPUTFILE}"
echo "  CHEM_CONFIG=${CHEM_CONFIG}"
echo "  OUTPUTFILE=${OUTPUTFILE}"
echo "  NUMTICKS=${NUMTICKS}  OMP_NUM_THREADS=${OMP_NUM_THREADS}  PROFILE=${IVDBM_PROFILE:-cpu}"

finalize_run_params() {
  local status="$1"
  local finished_at
  finished_at="$(date -Iseconds)"

  if [[ -f "${FINALIZE_SCRIPT}" && -f "${IVDBM_RUN_PARAMS_JSON}" ]]; then
    python3 "${FINALIZE_SCRIPT}" "${IVDBM_RUN_PARAMS_JSON}" \
      --timing "${RUN_TIMING_JSON}" \
      --time-file "${TIME_FILE}" \
      --sacct-job-id "${SLURM_JOB_ID:-}" \
      --started "${JOB_STARTED_AT}" \
      --finished "${finished_at}" \
      --hostname "${JOB_HOST}" \
      --profile "${IVDBM_PROFILE:-cpu}" \
      --exit-code "${status}" || \
      echo "testrun.sbatch: warning: could not update ${IVDBM_RUN_PARAMS_JSON}" >&2
  fi
}

set +e
/usr/bin/time -f 'wall_seconds=%e max_rss_kb=%M user_seconds=%U sys_seconds=%S' \
  -o "${TIME_FILE}" \
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
set -e

echo "testRun finished with exit code ${status}"
finalize_run_params "${status}"
exit "${status}"
