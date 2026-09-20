#!/bin/bash
# Slurm batch script for testRun — submit via scripts/submit_testrun.sh, not directly.
# Required env (set by submit_testrun.sh): REPO_ROOT, TESTRUN
# Optional env: OUTPUT_DIR, NUMTICKS, WXW, WYW, WZW, OUTPUTFILE, SIMULATION_CONFIG, IVDBM_PROFILE, KEEP_NAME, KEEP_WHY, RUN_COMMIT

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

mkdir -p "${OUTPUT_DIR}" logs || fail "cannot create ${OUTPUT_DIR} or logs/"

[[ -x "${TESTRUN:?Set TESTRUN when submitting this job.}" ]] \
  || fail "testRun not executable: ${TESTRUN}"
[[ -f "${SIMULATION_CONFIG}" ]] \
  || fail "simulation config missing: ${REPO_ROOT}/${SIMULATION_CONFIG} (copy configFiles/simulation_config.template.json)"

export LD_LIBRARY_PATH="${REPO_ROOT}/build/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
export IVDBM_OUTPUT_DIR="${OUTPUT_DIR}"
export IVDBM_RUN_PARAMS_JSON="${OUTPUT_DIR}/run_params.json"

if [[ -f /etc/profile.d/modules.sh ]]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
fi
module load StdEnv/2023 2>&1 || true
if [[ "${IVDBM_PROFILE:-cpu}" == "gpu" ]]; then
  module load cuda 2>&1 || module load cuda/12.2 2>&1 || true
fi

PARAVIEW_FLAG=()
if [[ "${PARAVIEW:-0}" -eq 1 ]]; then
  PARAVIEW_FLAG=(--paraview)
fi

echo "Job ${SLURM_JOB_ID:-?} starting on $(hostname) at $(date -Iseconds)"
echo "  REPO_ROOT=${REPO_ROOT}"
echo "  TESTRUN=${TESTRUN}"
echo "  OUTPUT_DIR=${OUTPUT_DIR}"
echo "  SIMULATION_CONFIG=${SIMULATION_CONFIG}"
echo "  OUTPUTFILE=${OUTPUTFILE}"
echo "  NUMTICKS=${NUMTICKS}  OMP_NUM_THREADS=${OMP_NUM_THREADS}  PROFILE=${IVDBM_PROFILE:-cpu}"

set +e
"${TESTRUN}" \
  --output-dir "${OUTPUT_DIR}" \
  --numticks "${NUMTICKS}" \
  --config "${SIMULATION_CONFIG}" \
  --outputfile "${OUTPUTFILE}" \
  --wxw "${WXW}" \
  --wyw "${WYW}" \
  --wzw "${WZW}" \
  "${PARAVIEW_FLAG[@]}"
status=$?
set -e

echo "testRun finished with exit code ${status}"

if [[ "${status}" -eq 0 && -n "${KEEP_NAME:-}" ]]; then
  RUN_FOLDER="${IVDBM_RESULTS_DIR}/$(date +%F)_${SLURM_JOB_ID}_${KEEP_NAME}"
  set +e
  mkdir -p "${RUN_FOLDER}/summary" \
    && cp -r "${OUTPUT_DIR}" "${RUN_FOLDER}/raw" \
    && cp "${OUTPUTFILE}" "${RUN_FOLDER}/summary/" \
    && cat > "${RUN_FOLDER}/manifest.txt" <<EOF
date: $(date '+%Y-%m-%d %H:%M')
cluster: ${CC_CLUSTER:-$(hostname)}
job_id: ${SLURM_JOB_ID}
commit: ${RUN_COMMIT:-unknown}
config: ${SIMULATION_CONFIG}
numticks: ${NUMTICKS}
world_mm: ${WXW}x${WYW}x${WZW}
why: ${KEEP_WHY:-}
EOF
  archive_status=$?
  set -e
  if [[ "${archive_status}" -eq 0 ]]; then
    echo "Archived run to ${RUN_FOLDER}"
  else
    echo "testrun.sbatch: warning: archiving to ${RUN_FOLDER} failed; output is still in ${OUTPUT_DIR}" >&2
  fi
fi

exit "${status}"
