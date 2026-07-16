# IVDBM-ABM

Agent-based model of stem cells in alginate hydrogel for intervertebral disc biomaterial studies. **Work in progress.**

## What this simulates

The model tracks stem cells growing inside an alginate hydrogel scaffold. The space is divided into small regions. Cells live in those regions, and the extracellular matrix (collagen, aggrecan, and related components) changes over time. Cytokines (TNF, TGF, IL-1β, O₂) spread between regions.

One tick is 30 minutes of simulated time. Set the tick length in the chemistry config file.

The main output is `Output_Biomarkers.csv`, written each tick with cell counts, cytokine levels, and scaffold measurements.

## First-time setup

You need CMake, a C++ compiler (GCC or Clang), and OpenGL/GLUT development libraries (on Debian/Ubuntu: `freeglut3-dev`). CUDA is optional and only needed for GPU diffusion builds.

From the repository root:

```bash
cp configFiles/chemical_environment.template.json configFiles/chemical_environment.json
mkdir -p output
```

Chemistry field details are in [configFiles/README.md](configFiles/README.md).

## What to change

| What you want to change | Where to set it | When |
| --- | --- | --- |
| Cytokine diffusivities, baselines, tick length | `configFiles/chemical_environment.json` | Before a run |
| Cell seeding, alginate, scaffold layout | `configFiles/config_scaffold.txt` | Before a run |
| Calibration coefficients | `Sample.txt` (repo root) | Before a run; not a command-line flag |
| Tick count, world size, paths, output | `testRun` flags or `submit_testrun.sh` flags | Each run |
| Shared DRAC simulation defaults | `scripts/job_defaults.env` | Optional |
| Model type, 3D, GPU diffusion, biomarker output | `src/common.h` | Rebuild after editing |
| Job time limit, CPUs, GPU choice | `submit_testrun.sh` flags | Each DRAC job |
| Unit test job defaults | `scripts/test_defaults.env` or `submit_tests.sh` | Test jobs only |

Command-line flags override `job_defaults.env`, which overrides the built-in defaults. With no flags on your own machine: 432 ticks, 3.1 mm cube. On DRAC with no extra flags: 24 ticks, 1 mm cube.

Biology settings (JSON, text configs, `Sample.txt`) are read when the run starts and are not affected by that order.

## Build

Binaries go to `build/bin/`. Libraries go to `build/lib/`.

Local machine:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SRC_TESTS=OFF
cmake --build build --target testRun -j"$(nproc)"
```

DRAC login node (load modules first; load CUDA only for GPU builds):

```bash
module load StdEnv/2023 gcc/12.3 cmake
./scripts/build_drac.sh

module load StdEnv/2023 gcc/12.3 cuda/12.2 cmake
./scripts/build_drac.sh --cuda --arch 90
./scripts/build_drac.sh --cuda --arch 90 --tests
```

For GPU builds you must load a CUDA module before running the script, or CMake will not find `nvcc`.

## Run

Edit the variables at the top of each block, then run from the repository root. One tick is 30 minutes of simulated time.

| Parameter | Flag or file | Notes |
| --- | --- | --- |
| Tick count | `--numticks` | 24 ≈ 12 h, 48 ≈ 1 day, 432 ≈ 9 days |
| World size (mm) | `--wxw`, `--wyw`, `--wzw` | Local default 3.1 mm; DRAC default 1 mm. GPU: 3.1 mm needs an H100 (~5 GB); use 1 mm or smaller on other GPUs |
| Cells / scaffold | `--inputfile` | e.g. `configFiles/config_scaffold.txt` |
| Chemistry | `--chem-config` | `configFiles/chemical_environment.json` |
| Calibration | `Sample.txt` | Repo root; not a command-line flag |
| Output directory | `--output-dir` | All CSVs; on DRAC also `run_params.json` |
| Biomarker CSV | `--outputfile` | Default `<output-dir>/Output_Biomarkers.csv` |

Local:

```bash
# --- simulation ---
NUMTICKS=48
WXW=3.1
WYW=3.1
WZW=3.1
INPUTFILE="configFiles/config_scaffold.txt"
CHEM_CONFIG="configFiles/chemical_environment.json"
OUTPUT_DIR="output/local_${NUMTICKS}ticks"
OUTPUTFILE="${OUTPUT_DIR}/Output_Biomarkers.csv"

./build/bin/testRun \
  --output-dir "${OUTPUT_DIR}" \
  --numticks "${NUMTICKS}" \
  --wxw "${WXW}" --wyw "${WYW}" --wzw "${WZW}" \
  --inputfile "${INPUTFILE}" \
  --chem-config "${CHEM_CONFIG}" \
  --outputfile "${OUTPUTFILE}"
```

DRAC (same simulation variables plus cluster options):

```bash
# --- simulation ---
EMAIL="$EMAIL"
NUMTICKS=48
WXW=3.1
WYW=3.1
WZW=3.1
INPUTFILE="configFiles/config_scaffold.txt"
CHEM_CONFIG="configFiles/chemical_environment.json"
OUTPUTFILE=""   # leave empty for default: <output-dir>/Output_Biomarkers.csv

# --- cluster ---
PROFILE="gpu"              # gpu after ./scripts/build_drac.sh --cuda; cpu for CPU-only build
WALLTIME="0-02:00:00"      # raise for long runs (submit default is 5 minutes)
CPUS=32
MEM="32000M"

SUBMIT=(./scripts/submit_testrun.sh "${EMAIL}" \
  --numticks "${NUMTICKS}" \
  --wxw "${WXW}" --wyw "${WYW}" --wzw "${WZW}" \
  --inputfile "${INPUTFILE}" \
  --chem-config "${CHEM_CONFIG}" \
  --profile "${PROFILE}" \
  --time "${WALLTIME}" \
  --cpus "${CPUS}" \
  --mem "${MEM}")

[[ -n "${OUTPUTFILE}" ]] && SUBMIT+=(--outputfile "${OUTPUTFILE}")

"${SUBMIT[@]}"
```

Append `--dry-run` to preview without submitting.

`./scripts/submit_testrun.sh "${EMAIL}" 12` is the same as `--numticks 12`.

To change shared DRAC simulation defaults, edit `scripts/job_defaults.env`. To change biology, edit files under `configFiles/` and `Sample.txt`.

## Results

| What | Local run | DRAC job |
| --- | --- | --- |
| Progress and errors | Terminal | `logs/testrun_<array>_<jobid>.out` and `.err` |
| CSVs and calibration data | `<output-dir>/` (default `output/`) | `output/job_<jobid>/` |
| Run settings snapshot | Only if `$IVDBM_RUN_PARAMS_JSON` is set | `<output-dir>/run_params.json` |

Main output files (default scaffold build):

| File | When | Contents |
| --- | --- | --- |
| `Output_Biomarkers.csv` | Each tick | Cell counts, cytokines, ECM, scaffold mechanics |
| `tgf_line.csv` | Each tick | TGF concentration along the x-axis |
| `SensitivityAnalysis/FinalTotalChemVR.dat` | End of run | Chemical totals for calibration |

After a DRAC job:

```bash
less logs/testrun_0_<jobid>.out
ls output/job_<jobid>/
cat output/job_<jobid>/run_params.json
```

## Testing on DRAC

Unit tests are separate from simulation runs. Build once with `--tests`, then submit test jobs as needed.

This project uses Catch2 v2. Catch2 v3 will not work. On some clusters the build script downloads Catch2 v2 during configure if no module is available.

```bash
./scripts/build_drac.sh --tests
./scripts/build_drac.sh --cuda --tests --arch 90
```

```bash
EMAIL="$EMAIL"

./scripts/submit_tests.sh "${EMAIL}" --suite cpu
./scripts/submit_tests.sh "${EMAIL}" --suite chemistry
./scripts/submit_tests.sh "${EMAIL}" --suite gpu --profile gpu
./scripts/submit_tests.sh "${EMAIL}" --suite all --time 1-00:00:00
./scripts/submit_tests.sh "${EMAIL}" --suite cpu --dry-run
```

| `--suite` | What it runs | GPU needed |
| --- | --- | --- |
| `chemistry` | Chemistry tests | No |
| `diffusion3d` | Diffusion3D tests | Mixed |
| `world` | World tests | No |
| `cpu` | All CPU-labelled tests | No |
| `gpu` | GPU-labelled tests | Yes |
| `all` | Full suite | Use CUDA build and `--profile gpu` for GPU coverage |

Test job logs: `logs/tests_<array>_<jobid>.out` / `.err`. CTest output: `output/tests_<jobid>/ctest.log`.

Shared test job defaults: `scripts/test_defaults.env`.

## Related docs

- [configFiles/README.md](configFiles/README.md): chemistry JSON and cell/scaffold text configs
- [src/Diffusion3D/README.md](src/Diffusion3D/README.md): diffusion module, tests, ParaView demo
