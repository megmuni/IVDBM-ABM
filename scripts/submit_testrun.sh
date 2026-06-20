#!/usr/bin/env bash
# Submit an IVDBM-ABM testRun job to a Digital Research Alliance of Canada (DRAC) cluster.
#
# Usage:
#   ./scripts/submit_testrun.sh user@institution.ca
#   ./scripts/submit_testrun.sh --email user@institution.ca --testrun build/bin/testRun

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SBATCH_SCRIPT="${SCRIPT_DIR}/testrun.sbatch"

EMAIL=""
ACCOUNT="def-nicoleli"
ARRAY="0-0"
TIME="0-00:05:00"
CPUS="32"
GPUS="h100:2"
MEM="32000M"
NUMTICKS="24"
WXW="1"
WYW="1"
WZW="1"
INPUTFILE="configFiles/config_scaffold.txt"
TESTRUN_REL="build/bin/testRun"
OUTPUT_DIR=""
DRY_RUN=0

usage() {
  cat <<'EOF'
Submit testRun to DRAC (Slurm).

Usage:
  submit_testrun.sh EMAIL [options]
  submit_testrun.sh --email EMAIL [options]

Required:
  EMAIL                 Address for Slurm mail notifications (--mail-user)

Options:
  --account ACCOUNT     Slurm allocation (default: def-nicoleli)
  --array RANGE         Job array range (default: 0-0, single task)
  --time D-HH:MM:SS     Wall time (default: 0-00:05:00)
  --cpus N              CPUs per task (default: 32)
  --gpus SPEC           GPUs per node (default: h100:2)
  --mem SIZE            Memory per node (default: 32000M)
  --testrun PATH        Path to testRun (default: build/bin/testRun, relative to repo root)
  --output-dir PATH     Override run output directory (default: output/job_<SLURM_JOB_ID> on compute node)
  --numticks N          Simulation ticks (default: 24, ~0.5 days at 30 min/tick)
  --wxw MM              World X width in mm (default: 1)
  --wyw MM              World Y width in mm (default: 1)
  --wzw MM              World Z width in mm (default: 1)
  --inputfile PATH      Input config relative to repo root
  --dry-run             Print sbatch command without submitting
  -h, --help            Show this help

The batch script is scripts/testrun.sbatch (do not sbatch it directly unless you
export REPO_ROOT, TESTRUN, etc. yourself). This wrapper sets those for you.

Examples:
  ./scripts/submit_testrun.sh user@institution.ca
  ./scripts/submit_testrun.sh user@institution.ca --testrun bin/testRun
  ./scripts/submit_testrun.sh user@institution.ca --array 0-0 --numticks 50
EOF
}

die() {
  echo "submit_testrun.sh: $*" >&2
  exit 1
}

valid_email() {
  [[ "$1" =~ ^[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}$ ]]
}

resolve_testrun() {
  local path="$1"
  if [[ "${path}" != /* ]]; then
    path="${REPO_ROOT}/${path}"
  fi
  [[ -f "${path}" && -x "${path}" ]] || die "testRun not found or not executable: ${path}"
  printf '%s' "${path}"
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
      shift 2
      ;;
    --mem)
      [[ $# -ge 2 ]] || die "missing value for $1"
      MEM="$2"
      shift 2
      ;;
    --testrun)
      [[ $# -ge 2 ]] || die "missing value for $1"
      TESTRUN_REL="$2"
      shift 2
      ;;
    --numticks)
      [[ $# -ge 2 ]] || die "missing value for $1"
      NUMTICKS="$2"
      shift 2
      ;;
    --wxw)
      [[ $# -ge 2 ]] || die "missing value for $1"
      WXW="$2"
      shift 2
      ;;
    --wyw)
      [[ $# -ge 2 ]] || die "missing value for $1"
      WYW="$2"
      shift 2
      ;;
    --wzw)
      [[ $# -ge 2 ]] || die "missing value for $1"
      WZW="$2"
      shift 2
      ;;
    --inputfile)
      [[ $# -ge 2 ]] || die "missing value for $1"
      INPUTFILE="$2"
      shift 2
      ;;
    --output-dir)
      [[ $# -ge 2 ]] || die "missing value for $1"
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
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

[[ -f "${SBATCH_SCRIPT}" ]] || die "missing batch script: ${SBATCH_SCRIPT}"

TESTRUN="$(resolve_testrun "${TESTRUN_REL}")"

[[ -f "${REPO_ROOT}/${INPUTFILE}" ]] || die "input file not found: ${REPO_ROOT}/${INPUTFILE}"

mkdir -p "${REPO_ROOT}/logs"

EXPORT_VARS="ALL,REPO_ROOT=${REPO_ROOT},TESTRUN=${TESTRUN},NUMTICKS=${NUMTICKS},WXW=${WXW},WYW=${WYW},WZW=${WZW},INPUTFILE=${INPUTFILE}"
if [[ -n "${OUTPUT_DIR}" ]]; then
  EXPORT_VARS="${EXPORT_VARS},OUTPUT_DIR=${OUTPUT_DIR}"
fi

SBATCH_ARGS=(
  --job-name=ivdbm-testrun
  --account="${ACCOUNT}"
  --mail-user="${EMAIL}"
  --mail-type=ALL
  --time="${TIME}"
  --nodes=1
  --cpus-per-task="${CPUS}"
  --gpus-per-node="${GPUS}"
  --mem="${MEM}"
  --array="${ARRAY}"
  --chdir="${REPO_ROOT}"
  --output="${REPO_ROOT}/logs/testrun_%a_%j.out"
  --error="${REPO_ROOT}/logs/testrun_%a_%j.err"
  --export="${EXPORT_VARS}"
  "${SBATCH_SCRIPT}"
)

if [[ "${DRY_RUN}" -eq 1 ]]; then
  printf 'Would run: sbatch'
  printf ' %q' "${SBATCH_ARGS[@]}"
  printf '\n'
  exit 0
fi

command -v sbatch >/dev/null 2>&1 || die "sbatch not found (log in to a DRAC cluster login node first)"

SUBMIT_OUTPUT="$(sbatch "${SBATCH_ARGS[@]}")"
echo "${SUBMIT_OUTPUT}"
JOB_ID="${SUBMIT_OUTPUT##* }"
echo "testRun: ${TESTRUN}"
if [[ -n "${OUTPUT_DIR}" ]]; then
  echo "Output dir (override): ${OUTPUT_DIR}"
else
  echo "Output dir: output/job_<SLURM_JOB_ID> (set on compute node)"
fi
echo "Mail notifications for ${EMAIL} will be sent by Slurm for job ${JOB_ID}."
