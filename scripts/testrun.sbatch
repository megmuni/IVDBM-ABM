#!/bin/bash
# Slurm batch script for testRun — submit via scripts/submit_testrun.sh, not directly.
# Required env (set by submit_testrun.sh): REPO_ROOT, TESTRUN
# Optional env: OUTPUT_DIR, NUMTICKS, WXW, WYW, WZW, INPUTFILE, OUTPUTFILE

set -euo pipefail

cd "${REPO_ROOT:?Set REPO_ROOT when submitting this job.}"

OUTPUT_DIR="${OUTPUT_DIR:-output/job_${SLURM_JOB_ID}}"
mkdir -p "${OUTPUT_DIR}" logs

module load StdEnv/2023 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12.2 2>/dev/null || true

INPUTFILE="${INPUTFILE:-configFiles/config_scaffold.txt}"
NUMTICKS="${NUMTICKS:-24}"
WXW="${WXW:-1}"
WYW="${WYW:-1}"
WZW="${WZW:-1}"
OUTPUTFILE="${OUTPUTFILE:-${OUTPUT_DIR}/Output_Biomarkers.csv}"

export IVDBM_OUTPUT_DIR="${OUTPUT_DIR}"

exec "${TESTRUN:?Set TESTRUN when submitting this job.}" \
  --output-dir "${OUTPUT_DIR}" \
  --numticks "${NUMTICKS}" \
  --inputfile "${INPUTFILE}" \
  --outputfile "${OUTPUTFILE}" \
  --wxw "${WXW}" \
  --wyw "${WYW}" \
  --wzw "${WZW}"
