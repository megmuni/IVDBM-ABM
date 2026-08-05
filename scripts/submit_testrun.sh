#!/usr/bin/env bash
# Submit an IVDBM-ABM testRun job to a Digital Research Alliance of Canada (DRAC) cluster.
#
# Usage:
#   ./scripts/submit_testrun.sh $EMAIL
#   ./scripts/submit_testrun.sh $EMAIL 48
#   ./scripts/submit_testrun.sh --email $EMAIL --profile gpu

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SBATCH_SCRIPT="${SCRIPT_DIR}/testrun.sbatch"
JOB_DEFAULTS="${SCRIPT_DIR}/job_defaults.env"

[[ -f "${JOB_DEFAULTS}" ]] || {
  echo "submit_testrun.sh: missing ${JOB_DEFAULTS}" >&2
  exit 1
}
# shellcheck source=job_defaults.env
source "${JOB_DEFAULTS}"

EMAIL=""
ACCOUNT="${IVDBM_DEFAULT_ACCOUNT}"
ARRAY="${IVDBM_DEFAULT_ARRAY}"
TIME="${IVDBM_DEFAULT_TIME}"
CPUS="${IVDBM_DEFAULT_CPUS}"
GPUS=""
MEM="${IVDBM_DEFAULT_MEM}"
PROFILE="${IVDBM_DEFAULT_PROFILE}"
NUMTICKS="${NUMTICKS}"
WXW="${WXW}"
WYW="${WYW}"
WZW="${WZW}"
SIMULATION_CONFIG="${SIMULATION_CONFIG}"
TESTRUN_REL="${IVDBM_DEFAULT_TESTRUN}"
OUTPUT_DIR=""
OUTPUTFILE=""
DRY_RUN=0
PARAVIEW=0
GPUS_EXPLICIT=0
TIME_EXPLICIT=0

usage() {
  cat <<'EOF'
Submit testRun to DRAC (Slurm).

Usage:
  submit_testrun.sh EMAIL [NUMTICKS] [options]
  submit_testrun.sh --email EMAIL [options]

Required:
  EMAIL                 Address for Slurm mail notifications (--mail-user)

Optional positional:
  NUMTICKS              Shorthand for --numticks (e.g. 48 for a 1-day run at 30 min/tick)

Options:
  --account ACCOUNT     Slurm allocation (default: def-nicoleli)
  --array RANGE         Job array range (default: 0-0, single task)
  --time D-HH:MM:SS     Wall time (default: 0-00:05:00)
  --cpus N              CPUs per task (default: 32)
  --profile cpu|gpu     Job profile (default: gpu — 2× H100; use cpu for CPU-only builds)
  --gpus SPEC           GPUs per node (overrides --profile; gpu default: h100:2)
  --mem SIZE            Memory per node (default: 32000M)
  --testrun PATH        Path to testRun (default: build/bin/testRun, relative to repo root)
  --output-dir PATH     Override run output directory (default: output/job_<SLURM_JOB_ID> on compute node)
  --outputfile PATH     Biomarker CSV path (default: <output-dir>/Output_Biomarkers.csv on compute node)
  --config PATH         Simulation JSON (default: configFiles/simulation_config.json)
  --numticks N          Simulation ticks (default: 24, ~12 h at 30 min/tick)
  --wxw MM              World X width in mm (default: 1)
  --wyw MM              World Y width in mm (default: 1)
  --wzw MM              World Z width in mm (default: 1)
  --paraview            Writes <output-dir>/paraview/patches_t*.vti, chem_t*.vti, and ecm_t*.vti every 12 ticks
  --dry-run             Print sbatch command without submitting
  -h, --help            Show this help

Defaults live in scripts/job_defaults.env. The batch script is scripts/testrun.sbatch
(submit via this wrapper, not sbatch directly unless you export REPO_ROOT, TESTRUN, etc.).

Examples:
  ./scripts/submit_testrun.sh $EMAIL
  ./scripts/submit_testrun.sh $EMAIL 48
  ./scripts/submit_testrun.sh $EMAIL --profile gpu --numticks 24
  ./scripts/submit_testrun.sh $EMAIL --config configFiles/simulation_config.json
EOF
}

die() {
  echo "submit_testrun.sh: $*" >&2
  exit 1
}

warn() {
  echo "submit_testrun.sh: warning: $*" >&2
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

testrun_links_cuda() {
  ldd "${TESTRUN}" 2>/dev/null | grep -qE 'libcudart|libcuda'
}

warn_profile_binary_mismatch() {
  if testrun_links_cuda; then
    if [[ "${PROFILE}" == "cpu" ]]; then
      warn "testRun is CUDA-linked; default is --profile gpu (pass --profile gpu explicitly)"
    fi
  elif [[ "${PROFILE}" == "gpu" ]]; then
    warn "profile=gpu but testRun is not CUDA-linked; run ./scripts/build_drac.sh --cuda or use --profile cpu"
  fi
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
      TIME_EXPLICIT=1
      shift 2
      ;;
    --cpus)
      [[ $# -ge 2 ]] || die "missing value for $1"
      CPUS="$2"
      shift 2
      ;;
    --profile)
      [[ $# -ge 2 ]] || die "missing value for $1"
      PROFILE="$2"
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
    --config)
      [[ $# -ge 2 ]] || die "--config requires a path"
      SIMULATION_CONFIG="$2"
      shift 2
      ;;
    --outputfile)
      [[ $# -ge 2 ]] || die "missing value for $1"
      OUTPUTFILE="$2"
      shift 2
      ;;
    --output-dir)
      [[ $# -ge 2 ]] || die "missing value for $1"
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --paraview)
      PARAVIEW=1
      shift
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
        if [[ "$1" =~ ^[0-9]+$ ]]; then
          NUMTICKS="$1"
          shift
        else
          die "unexpected argument: $1 (email already set to ${EMAIL}; use --numticks N or pass N as second argument)"
        fi
      else
        EMAIL="$1"
        shift
      fi
      ;;
  esac
done

[[ -n "${EMAIL}" ]] || die "EMAIL is required (positional or --email)"
valid_email "${EMAIL}" || die "invalid email address: ${EMAIL}"

if [[ "${GPUS_EXPLICIT}" -eq 0 ]]; then
  apply_profile
fi

if [[ "${NUMTICKS}" -gt 24 && "${TIME_EXPLICIT}" -eq 0 && "${TIME}" == "${IVDBM_DEFAULT_TIME}" ]]; then
  warn "NUMTICKS=${NUMTICKS} with default wall time ${TIME} may TIMEOUT; pass --time D-HH:MM:SS"
fi

[[ -f "${SBATCH_SCRIPT}" ]] || die "missing batch script: ${SBATCH_SCRIPT}"

TESTRUN="$(resolve_testrun "${TESTRUN_REL}")"

if [[ ! -f "${REPO_ROOT}/${SIMULATION_CONFIG}" ]]; then
  die "missing ${SIMULATION_CONFIG} — run: cp configFiles/simulation_config.template.json configFiles/simulation_config.json"
fi

warn_profile_binary_mismatch

mkdir -p "${REPO_ROOT}/logs"

EXPORT_VARS="NONE,REPO_ROOT=${REPO_ROOT},TESTRUN=${TESTRUN},NUMTICKS=${NUMTICKS},WXW=${WXW},WYW=${WYW},WZW=${WZW},SIMULATION_CONFIG=${SIMULATION_CONFIG},IVDBM_PROFILE=${PROFILE}"
if [[ -n "${OUTPUT_DIR}" ]]; then
  EXPORT_VARS="${EXPORT_VARS},OUTPUT_DIR=${OUTPUT_DIR}"
fi
if [[ -n "${OUTPUTFILE}" ]]; then
  EXPORT_VARS="${EXPORT_VARS},OUTPUTFILE=${OUTPUTFILE}"
fi

SBATCH_ARGS=(
  --job-name=ivdbm-testrun
  --account="${ACCOUNT}"
  --mail-user="${EMAIL}"
  --mail-type=ALL
  --time="${TIME}"
  --nodes=1
  --cpus-per-task="${CPUS}"
  --mem="${MEM}"
  --array="${ARRAY}"
  --chdir="${REPO_ROOT}"
  --output="${REPO_ROOT}/logs/testrun_%a_%j.out"
  --error="${REPO_ROOT}/logs/testrun_%a_%j.err"
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
  exit 0
fi

command -v sbatch >/dev/null 2>&1 || die "sbatch not found (log in to a DRAC cluster login node first)"

SUBMIT_OUTPUT="$(sbatch "${SBATCH_ARGS[@]}")"
echo "${SUBMIT_OUTPUT}"
JOB_ID="${SUBMIT_OUTPUT##* }"
echo "testRun: ${TESTRUN}"
echo "Profile: ${PROFILE}${GPUS:+ (${GPUS})}"
if [[ -n "${OUTPUT_DIR}" ]]; then
  echo "Output dir (override): ${OUTPUT_DIR}"
else
  echo "Output dir: output/job_<SLURM_JOB_ID> (set on compute node)"
fi
echo "Mail notifications for ${EMAIL} will be sent by Slurm for job ${JOB_ID}."
