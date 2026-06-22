# IVDBM-ABM

Agent-based model (C++) of stem cell–seeded alginate hydrogel for IVD biomaterial studies. **Work in progress.**

## Prerequisites

- CMake 3.x, a C++ compiler (GCC or Clang; C++11 minimum, C++17 for Chemistry/Diffusion3D)
- OpenGL / GLU / GLUT development libraries (e.g. on Debian/Ubuntu: `freeglut3-dev`)
- Optional: CUDA toolkit and GPU nodes when building with GPU diffusion (`-DDIFFUSION3D_CUDA=ON`)

Before the first run, create the chemical environment config (see `[configFiles/README.md](configFiles/README.md)`):

```bash
cp configFiles/chemical_environment.template.json configFiles/chemical_environment.json
mkdir -p output
```

Simulation and cluster parameters are documented in [Where to set run parameters](#where-to-set-run-parameters) (including [`scripts/job_defaults.env`](scripts/job_defaults.env) for shared DRAC defaults).

## Where to set run parameters

Run parameters live in **different places** depending on what you are changing and how you launch the model. There is no single `.env` file for the whole project.

### What to edit when

| You want to change… | Set it here | When |
| ------------------- | ----------- | ---- |
| Tick count, world size, config paths, output location for **one run** | `testRun` CLI flags (local) or [`scripts/submit_testrun.sh`](scripts/submit_testrun.sh) flags (DRAC) | Every one-off or experimental run |
| Shared **simulation** defaults for all DRAC jobs in your group | [`scripts/job_defaults.env`](scripts/job_defaults.env) (`NUMTICKS`, `WXW`, paths, …) | When the whole group should stop passing the same flags |
| Slurm wall time, CPUs, memory, account, GPU profile | `submit_testrun.sh` flags (`--time`, `--cpus`, …) or `IVDBM_DEFAULT_*` in `job_defaults.env` | Per job, or group-wide cluster policy |
| Chemistry (diffusivities, baselines, tick interval) | [`configFiles/chemical_environment.json`](configFiles/chemical_environment.json) | Before any run; copy from [`configFiles/chemical_environment.template.json`](configFiles/chemical_environment.template.json) once |
| Cell seeding, alginate, scaffold layout | [`configFiles/config_scaffold.txt`](configFiles/config_scaffold.txt) (or `--inputfile`) | Before any run |
| Calibration coefficients | [`Sample.txt`](Sample.txt) at repo root | Before any run; **not** a CLI flag |
| Model type (scaffold vs vocal fold), 3D, GPU diffusion, biomarker output | [`src/common.h`](src/common.h) | **Recompile** after editing |
| CTest suite / test-job Slurm defaults | [`scripts/test_defaults.env`](scripts/test_defaults.env) or [`scripts/submit_tests.sh`](scripts/submit_tests.sh) | Test jobs only; unrelated to `testRun` simulation |

### Precedence (same knob, multiple sources)

For ticks, world size, and paths, **highest priority wins**:

```text
submit_testrun.sh / testRun CLI  >  job_defaults.env  >  C++ defaults (input_utils.h)
```

- **Local:** `./build/bin/testRun` with no flags uses C++ defaults (**432 ticks**, **3.1 mm** cube).
- **DRAC:** [`scripts/job_defaults.env`](scripts/job_defaults.env) shortens that to **24 ticks**, **1 mm** cube for smoke tests unless you pass flags or edit the file.
- **Biology** (JSON, text configs, `Sample.txt`) is read at runtime and is independent of the precedence chain above.

### Default values (reference)

| Setting | Local (`testRun`, no flags) | DRAC (`submit_testrun.sh`, no extra flags) |
| ------- | ----------------------------- | -------------------------------------------- |
| `--numticks` | **432** (~9 days) | **24** (~12 h) via `job_defaults.env` |
| `--wxw`, `--wyw`, `--wzw` | **3.1 mm** each | **1 mm** each |
| `--inputfile` | `configFiles/config_scaffold.txt` | same |
| `--chem-config` | `configFiles/chemical_environment.json` | same |
| `--output-dir` | `output` (or `$IVDBM_OUTPUT_DIR`) | `output/job_<jobid>` on compute node |
| `--outputfile` | `<output-dir>/Output_Biomarkers.csv` | same |
| Tick duration | **30 min** (`tick_interval_minutes` in JSON) | same |
| `--patchwidth` | **0.01 mm** (fixed in code; not on CLI) | same |
| Calibration | **`Sample.txt`** (auto-loaded) | same |

Compile-time flags in [`src/common.h`](src/common.h) for the default scaffold build include `MODEL_SCAFFOLD`, `MODEL_3D`, `PDE_DIFFUSE`, `PEPTIDE_BM`, `BIOMARKER_OUTPUT`, and `CALIBRATION`.

### Copy-paste parameter blocks

Edit the variables at the top of each block, then run from the **repository root**. One tick = **30 minutes** simulated time.

| Parameter | Flag / file | Notes |
| --------- | ----------- | ----- |
| Tick count | `--numticks` | 24 ≈ 12 h, 48 ≈ 1 day, 432 ≈ 9 days |
| World size (mm) | `--wxw`, `--wyw`, `--wzw` | C++ default 3.1 mm; DRAC shared default 1 mm |
| Cells / scaffold | `--inputfile` | e.g. `configFiles/config_scaffold.txt` |
| Chemistry | `--chem-config` | `configFiles/chemical_environment.json` |
| Calibration | `Sample.txt` | Repo root; not a CLI flag |
| Output directory | `--output-dir` | All CSVs; on DRAC also `run_params.json` |
| Biomarker CSV | `--outputfile` | Default `<output-dir>/Output_Biomarkers.csv` |

**Local run:**

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

**DRAC (Slurm) - same simulation vars plus cluster options:**

```bash
# --- simulation ---
EMAIL="your.email@mail.mcgill.ca"
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

Preview without submitting: append `--dry-run` to the submit command.

**Shorthand (ticks only):** `./scripts/submit_testrun.sh "${EMAIL}" 12` is the same as `--numticks 12`.

**Prefer editing a file over flags?** Change simulation defaults in [`scripts/job_defaults.env`](scripts/job_defaults.env). Change biology in `configFiles/` and `Sample.txt` (see table above).

## Compile

Configure and build from the repository root. Binaries are written to `build/bin/` and libraries to `build/lib/`.

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SRC_TESTS=OFF
cmake --build build --target testRun -j"$(nproc)"
```

### DRAC (Alliance Canada)

On a cluster login node, load compiler modules before `cmake`. **If you pass `-DDIFFUSION3D_CUDA=ON`, you must load CUDA first**, otherwise CMake fails with `Failed to find nvcc`.

```bash
# CPU-only (default; no nvcc required)
module load StdEnv/2023 gcc/12.3 cmake
./scripts/build_drac.sh

# GPU diffusion (H100 = sm_90) + tests
module load StdEnv/2023 gcc/12.3 cuda/12.2 cmake
./scripts/build_drac.sh --cuda --arch 90 --tests
```

Or manually:

```bash
module load StdEnv/2023 gcc/12.3 cuda/12.2 cmake   # cuda only needed for GPU build
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SRC_TESTS=OFF \
  -DDIFFUSION3D_CUDA=ON \
  -DDIFFUSION3D_CUDA_ARCHITECTURES=90
cmake --build build --target testRun -j"$(nproc)"
```

To enable GPU diffusion kernels locally (requires CUDA toolkit / `nvcc` on `PATH`):

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SRC_TESTS=OFF \
  -DDIFFUSION3D_CUDA=ON \
  -DDIFFUSION3D_CUDA_ARCHITECTURES=90   # H100; use 89 for Ada (RTX 4050), etc.
cmake --build build --target testRun -j"$(nproc)"
```

Legacy CMake entry points (`CMakeListsSerial.txt`, `CMakeListsOMP.txt`, and Intel variants) remain for older cluster setups that linked the deprecated `Diffusion` library directly against CUDA. Prefer the root `CMakeLists.txt` above for current development.

Simulation builds keep tests off by default (`BUILD_SRC_TESTS=OFF`). Use `./scripts/build_drac.sh --tests` on DRAC to build test binaries (see [Testing on DRAC](#testing-on-drac)).

## Testing on DRAC

Catch2/CTest suites live under [`src/tests/`](src/tests/). They are **separate from `testRun` simulation jobs** - build once, then submit test jobs repeatedly via Slurm.

### Build test binaries

This project uses **Catch2 v2 only** (v3 is incompatible). On many DRAC clusters a `catch2/2.x` module is available; on others (e.g. Nibi) `./scripts/build_drac.sh --tests` falls back to downloading Catch2 v2.13.10 during CMake configure. Do not `module load catch2/3.x`.

```bash
# CPU tests (login node)
./scripts/build_drac.sh --tests

# CPU + GPU tests (requires cuda module)
./scripts/build_drac.sh --cuda --tests --arch 90
```

If the build warns that no Catch2 module was found, that is expected on some clusters — the FetchContent download runs once at configure time on the login node.

### Submit test jobs (Slurm)

```bash
EMAIL="your.email@mail.mcgill.ca"

# Fast CPU suite (default) - chemistry, diffusion3d CPU, world
./scripts/submit_tests.sh "${EMAIL}" --suite cpu

# Single module
./scripts/submit_tests.sh "${EMAIL}" --suite chemistry

# GPU-labelled tests (CUDA build + GPU node)
./scripts/submit_tests.sh "${EMAIL}" --suite gpu --profile gpu

# Everything
./scripts/submit_tests.sh "${EMAIL}" --suite all --time 1-00:00:00

# Preview
./scripts/submit_tests.sh "${EMAIL}" --suite cpu --dry-run
```

| `--suite`     | CTest filter     | GPU needed                                              |
| ------------- | ---------------- | ------------------------------------------------------- |
| `chemistry`   | `-L chemistry`   | No                                                      |
| `diffusion3d` | `-L diffusion3d` | Mixed (CPU tests always; GPU test if CUDA build)        |
| `world`       | `-L world`       | No                                                      |
| `cpu`         | `-L cpu`         | No                                                      |
| `gpu`         | `-L gpu`         | Yes                                                     |
| `all`         | (no filter)      | For full GPU coverage, use CUDA build + `--profile gpu` |

Pass extra ctest flags after `--`:

```bash
./scripts/submit_tests.sh "${EMAIL}" --suite chemistry -- --verbose
```

### Test output

| Artifact     | Location                                  |
| ------------ | ----------------------------------------- |
| Job log      | `logs/tests_<array>_<jobid>.out` / `.err` |
| JUnit XML    | `output/tests_<jobid>/test_results.xml`   |
| Summary JSON | `output/tests_<jobid>/test_summary.json`  |

Shared Slurm defaults: [`scripts/test_defaults.env`](scripts/test_defaults.env).

## Run `testRun` locally

Run from the **repository root** so default paths such as `configFiles/chemical_environment.json` resolve correctly. For where each parameter is set and copy-paste blocks, see [Where to set run parameters](#where-to-set-run-parameters).

```bash
./build/bin/testRun \
  --output-dir output/local_run \
  --numticks 24 \
  --inputfile configFiles/config_scaffold.txt \
  --wxw 1 --wyw 1 --wzw 1
```

If `--outputfile` is omitted, the biomarker CSV defaults to `<output-dir>/Output_Biomarkers.csv`.

Locally, `testRun` prints progress to the terminal; simulation files go under `--output-dir`. See [Logging and output](#logging-and-output) for a full list of files and when `run_params.json` is written.

Useful CLI flags (see `./build/bin/testRun --help`):

| Flag                      | Default                                 | Description                                              |
| ------------------------- | --------------------------------------- | -------------------------------------------------------- |
| `--numticks`              | 432                                     | Number of simulation ticks (1 tick = 30 min)             |
| `--wxw`, `--wyw`, `--wzw` | 3.1 mm                                  | World width, length, height                              |
| `--inputfile`             | `configFiles/config_scaffold.txt`       | Text config for cells, alginate, and scaffold parameters |
| `--outputfile`            | `<output-dir>/Output_Biomarkers.csv`    | CSV biomarker output path                                |
| `--output-dir`            | `output`                                | Base directory for all run artifacts                     |
| `--chem-config`           | `configFiles/chemical_environment.json` | Path to chemical environment JSON                        |

With `MODEL_SCAFFOLD` defined in `src/common.h`, defaults match a small alginate scaffold (`configFiles/config_scaffold.txt`, 3.1 mm cube unless overridden on the command line). Calibration values are read from **`Sample.txt`** in the repo root (not configurable via CLI). See [Where to set run parameters](#where-to-set-run-parameters).

## Logging and output

Every run produces information in up to **three places**. What you get depends on whether you run locally or submit to Slurm.

### Where things go

| Channel               | Local run                                           | Slurm job (DRAC)                                     |
| --------------------- | --------------------------------------------------- | ---------------------------------------------------- |
| **Console / job log** | Terminal stdout/stderr                              | `logs/testrun_<array>_<jobid>.out` and `.err`        |
| **Simulation files**  | `<output-dir>/` (default `output/`)                 | `output/job_<jobid>/` (override with `--output-dir`) |
| **Run manifest**      | Not written unless you set `$IVDBM_RUN_PARAMS_JSON` | `<output-dir>/run_params.json`                       |

Set the simulation output directory with `--output-dir` or `$IVDBM_OUTPUT_DIR`. On Slurm, the submit script passes this through automatically.

### Slurm job logs (`logs/`)

Submitted via [`scripts/submit_testrun.sh`](scripts/submit_testrun.sh):

- **`logs/testrun_<array>_<jobid>.out`** - combined stdout from the batch script and `testRun`
- **`logs/testrun_<array>_<jobid>.err`** - stderr (errors, some C++ warnings)

The `.out` log typically contains:

1. **Job startup** - host, paths, tick count, profile, `OMP_NUM_THREADS`
2. **`testRun` stdout** - resolved CLI options, world setup, hydrogel properties, per-tick agent/ECM counts (when `BIOMARKER_OUTPUT` is enabled), calibration parameter dumps (when `CALIBRATION` is enabled)
3. **Job finish** - `testRun` exit code and a **Job summary** block (wall time, max memory, average tick time, Slurm accounting when available)

Slurm also sends email notifications to the address passed to the submit script (`--mail-type=ALL`).

### Simulation output files (`<output-dir>/`)

Written under `--output-dir` (or `$IVDBM_OUTPUT_DIR`). With the default scaffold build (`BIOMARKER_OUTPUT` and `CALIBRATION` in [`src/common.h`](src/common.h)):

| File                                           | When       | Contents                                                                                                                            |
| ---------------------------------------------- | ---------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| **`Output_Biomarkers.csv`**                    | Each tick  | Clock, simulated day, cytokine totals, ECM, cell counts, scaffold mechanics (E, swelling, mass loss, viability, differentiation, …) |
| **`tgf_line.csv`**                             | Each tick  | TGF concentration along the x-axis (diffusion debugging)                                                                            |
| **`SensitivityAnalysis/FinalTotalChemVR.dat`** | End of run | Total chemical volumes/ratios for calibration workflows                                                                             |

Override the biomarker path with `--outputfile` (submit script or `testRun` CLI).

Optional compile-time outputs (off by default): ParaView VTK dumps under `Simulation/` when `PARAVIEW_RENDERING` is defined.

### Run manifest (`run_params.json`)

**Slurm jobs only** (unless you export `IVDBM_RUN_PARAMS_JSON` yourself). Written to `<output-dir>/run_params.json`.

| Phase                                                                            | Sections written                                                                    |
| -------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| **Job start** (`testRun` startup)                                                | `slurm` - job id, partition, CPUs, memory, account, profile, …                      |
|                                                                                  | `simulation` - resolved ticks, world size, config paths, calibration file           |
|                                                                                  | `command` - full `argv`                                                             |
| **Job end** ([`scripts/finalize_run_params.py`](scripts/finalize_run_params.py)) | `runtime` - start/finish timestamps, hostname, exit code, profile                   |
|                                                                                  | `runtime.testrun` - wall time and max RSS from `/usr/bin/time`                      |
|                                                                                  | `runtime.tick_timing` - setup time and per-tick execution stats from C++            |
|                                                                                  | `runtime.slurm_accounting` - elapsed time, CPU time, max memory, state from `sacct` |

After a job finishes, check **`run_params.json`** for a machine-readable record of what ran and how long it took. The same summary is printed at the bottom of the `.out` log.

### Quick reference after a DRAC job

```bash
# Job log (stdout)
less logs/testrun_0_<jobid>.out

# Results CSVs and manifest
ls output/job_<jobid>/
cat output/job_<jobid>/run_params.json
```

Override the output directory on submit: `./scripts/submit_testrun.sh email@host --output-dir output/my_run_<jobid>`

## Run `testRun` on DRAC (Slurm)

Use [`scripts/submit_testrun.sh`](scripts/submit_testrun.sh) from the repository root. It wraps [`scripts/testrun.sbatch`](scripts/testrun.sbatch) with Slurm options and paths. For parameter locations and a copy-paste block, see [Where to set run parameters](#where-to-set-run-parameters).

```bash
./scripts/submit_testrun.sh your.email@mail.mcgill.ca
```

Common options:

| Flag                 | Default                                 | Purpose                                        |
| -------------------- | --------------------------------------- | ---------------------------------------------- |
| `--testrun PATH`     | `build/bin/testRun`                     | Executable (relative to repo root or absolute) |
| `--profile cpu\|gpu` | `gpu` (2× H100)                         | Use `cpu` only for CPU-only builds             |
| `--array`            | `0-0`                                   | Single job (use e.g. `0-4` for a batch)        |
| `--account`          | `def-nicoleli`                          | Slurm allocation                               |
| `--numticks`         | `24` (~12 h simulated)                  | Simulation length (30 min/tick)                |
| `--chem-config`      | `configFiles/chemical_environment.json` | Chemical environment JSON                      |
| `--outputfile`       | `<output-dir>/Output_Biomarkers.csv`    | Biomarker CSV path                             |
| `--output-dir`       | `output/job_<jobid>`                    | Simulation output directory                    |
| `--inputfile`        | `configFiles/config_scaffold.txt`       | Cell/scaffold config                           |

Shared cluster defaults are in [`scripts/job_defaults.env`](scripts/job_defaults.env) (simulation + Slurm). Default submit profile is **gpu** (matches `./scripts/build_drac.sh --cuda`). Use `--profile cpu` for CPU-only builds. See [Where to set run parameters](#where-to-set-run-parameters) for when to edit that file vs passing flags.

Examples:

```bash
# Default binary at build/bin/testRun
./scripts/submit_testrun.sh your.email@mail.mcgill.ca

# In-tree binary from an older layout
./scripts/submit_testrun.sh your.email@mail.mcgill.ca --testrun bin/testRun

# Preview the sbatch command
./scripts/submit_testrun.sh your.email@mail.mcgill.ca --dry-run
```

See [Logging and output](#logging-and-output) for job logs, CSV artifacts, and `run_params.json`.

## Related docs

- `[configFiles/README.md](configFiles/README.md)` - chemical environment JSON and text config files
- `[src/Diffusion3D/README.md](src/Diffusion3D/README.md)` - Diffusion3D module, tests, and ParaView demo
