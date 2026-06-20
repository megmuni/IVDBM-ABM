#!/bin/bash
#SBATCH --time=0-00:05:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=h100:2
#SBATCH --mem=32000M
# Account, mail-user, array range, and log path are set by scripts/submit_testrun.sh.

set -euo pipefail

cd "${REPO_ROOT:?Set REPO_ROOT when submitting this job.}"

mkdir -p output logs

module load StdEnv/2023 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12.2 2>/dev/null || true

INPUTFILE="${INPUTFILE:-configFiles/config_scaffold.txt}"
NUMTICKS="${NUMTICKS:-288}"
WXW="${WXW:-1}"
WYW="${WYW:-1}"
WZW="${WZW:-1}"
OUTPUTFILE="${OUTPUTFILE:-output/Output_Biomarkers_${SLURM_ARRAY_TASK_ID}.csv}"

exec "${TESTRUN:?Set TESTRUN when submitting this job.}" \
  --numticks "${NUMTICKS}" \
  --inputfile "${INPUTFILE}" \
  --outputfile "${OUTPUTFILE}" \
  --wxw "${WXW}" \
  --wyw "${WYW}" \
  --wzw "${WZW}"
