# IVDBM-ABM

Agent-based model (C++) of stem cell–seeded alginate hydrogel for IVD biomaterial studies. **Work in progress.**

## Prerequisites

- CMake 3.x, a C++ compiler (GCC or Clang; C++11 minimum, C++17 for Chemistry/Diffusion3D)
- OpenGL / GLU / GLUT development libraries (e.g. on Debian/Ubuntu: `freeglut3-dev`)
- Optional: CUDA toolkit and GPU nodes when building with GPU diffusion (`-DDIFFUSION3D_CUDA=ON`)

Before the first run, create the chemical environment config (see [`configFiles/README.md`](configFiles/README.md)):

```bash
cp configFiles/chemical_environment.template.json configFiles/chemical_environment.json
mkdir -p output
```

Model compile-time options (scaffold vs vocal fold, 3D, biomarker output, etc.) are set in [`src/common.h`](src/common.h).

## Default run definition

These are the **canonical defaults** when you run `./build/bin/testRun` with no CLI flags from the repository root (with the current `src/common.h` scaffold build). Slurm submission via [`scripts/submit_testrun.sh`](scripts/submit_testrun.sh) passes shorter overrides (24 ticks, 1 mm cube) suited to quick cluster smoke tests; override those flags on the submit script if you want a full-length run.

| Setting                   | Default                                                                 |
| ------------------------- | ----------------------------------------------------------------------- |
| `--numticks`              | **432** (~9 days at 30 min/tick)                                        |
| `--wxw`, `--wyw`, `--wzw` | **3.1 mm** each (3D scaffold cube)                                      |
| `--patchwidth`            | **0.01 mm** (fixed; not overridable on CLI)                             |
| `--inputfile`             | `configFiles/config_scaffold.txt`                                       |
| `--chem-config`           | `configFiles/chemical_environment.json` (copy from template)            |
| `--output-dir`            | `output` (or `$IVDBM_OUTPUT_DIR`)                                       |
| `--outputfile`            | `<output-dir>/Output_Biomarkers.csv`                                    |
| Tick duration             | **30 minutes** per tick (ABM clock and `tick_interval_minutes` in JSON) |
| Calibration parameters    | **`Sample.txt`** at repo root (loaded automatically; not a CLI flag)    |

Compile-time flags in [`src/common.h`](src/common.h) for the default scaffold build include: `MODEL_SCAFFOLD`, `MODEL_3D`, `PDE_DIFFUSE`, `PEPTIDE_BM`, `BIOMARKER_OUTPUT`, and `CALIBRATION`. Change these before rebuilding to switch model geometry, diffusion backend, or output behaviour.

Biological chemistry (diffusivities, baselines, tick interval) lives in `configFiles/chemical_environment.json`; the template sets `tick_interval_minutes` to **30**.

### Adjusting simulation length

One tick = **30 minutes** of simulated time. Set `--numticks N` on `testRun` locally, or pass `--numticks N` (or a second positional argument) to [`scripts/submit_testrun.sh`](scripts/submit_testrun.sh) on DRAC.

| `--numticks` | Simulated time        |
| ------------ | --------------------- |
| 24           | 12 hours              |
| 48           | 1 day                 |
| 432          | 9 days (full default) |

**Local** (from repo root):

```bash
# Quick test — 24 ticks (~12 h simulated), small 1 mm cube
./build/bin/testRun \
  --output-dir output/local_24 \
  --numticks 24 \
  --wxw 1 --wyw 1 --wzw 1

# Full default — 432 ticks (~9 days simulated), 3.1 mm scaffold cube
./build/bin/testRun \
  --output-dir output/local_432 \
  --numticks 432
```

**DRAC (Slurm)**:

```bash
# Shorthand: email then tick count (same as --numticks 12)
./scripts/submit_testrun.sh your.email@mail.mcgill.ca 12

# GPU diffusion (build first: ./scripts/build_drac.sh --cuda)
./scripts/submit_testrun.sh your.email@mail.mcgill.ca --profile gpu --numticks 24

# Long run — also raise wall time (submit default is 5 minutes)
./scripts/submit_testrun.sh your.email@mail.mcgill.ca \
  --numticks 432 \
  --time 2-00:00:00
```

Cluster simulation defaults (24 ticks, 1 mm cube, etc.) live in [`scripts/job_defaults.env`](scripts/job_defaults.env). Edit that file or override flags on the submit script.

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

# GPU diffusion (H100 = sm_90)
module load StdEnv/2023 gcc/12.3 cuda/12.2 cmake
./scripts/build_drac.sh --cuda --arch 90
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

## Run `testRun` locally

Run from the **repository root** so default paths such as `configFiles/chemical_environment.json` resolve correctly.

```bash
./build/bin/testRun \
  --output-dir output/local_run \
  --numticks 24 \
  --inputfile configFiles/config_scaffold.txt \
  --wxw 1 --wyw 1 --wzw 1
```

If `--outputfile` is omitted, the biomarker CSV defaults to `<output-dir>/Output_Biomarkers.csv`.

Useful CLI flags (see `./build/bin/testRun --help`):

| Flag                      | Default                                 | Description                                              |
| ------------------------- | --------------------------------------- | -------------------------------------------------------- |
| `--numticks`              | 432                                     | Number of simulation ticks (1 tick = 30 min)             |
| `--wxw`, `--wyw`, `--wzw` | 3.1 mm                                  | World width, length, height                              |
| `--inputfile`             | `configFiles/config_scaffold.txt`       | Text config for cells, alginate, and scaffold parameters |
| `--outputfile`            | `<output-dir>/Output_Biomarkers.csv`    | CSV biomarker output path                                |
| `--output-dir`            | `output`                                | Base directory for all run artifacts                     |
| `--chem-config`           | `configFiles/chemical_environment.json` | Path to chemical environment JSON                        |

With `MODEL_SCAFFOLD` defined in `src/common.h`, defaults match a small alginate scaffold (`configFiles/config_scaffold.txt`, 3.1 mm cube unless overridden on the command line). Calibration values are read from **`Sample.txt`** in the repo root (not configurable via CLI). See [Default run definition](#default-run-definition).

## Run `testRun` on DRAC (Slurm)

Use [`scripts/submit_testrun.sh`](scripts/submit_testrun.sh) from the repository root. It calls [`scripts/testrun.sbatch`](scripts/testrun.sbatch) with the correct `#SBATCH` options, your email, and paths to the binary and config files.

```bash
./scripts/submit_testrun.sh your.email@mail.mcgill.ca
```

Common options:

| Flag               | Default                           | Purpose                                        |
| ------------------ | --------------------------------- | ---------------------------------------------- |
| `--testrun PATH`   | `build/bin/testRun`               | Executable (relative to repo root or absolute) |
| `--profile cpu\|gpu` | `gpu` (2× H100)                   | Use `cpu` only for `./scripts/build_drac.sh` without `--cuda` |
| `--array`          | `0-0`                             | Single job (use e.g. `0-4` for a batch)        |
| `--account`        | `def-nicoleli`                    | Slurm allocation                               |
| `--numticks`       | `24` (~12 h simulated)            | Simulation length (30 min/tick)                |
| `--chem-config`    | `configFiles/chemical_environment.json` | Chemical environment JSON              |
| `--outputfile`     | `<output-dir>/Output_Biomarkers.csv` | Biomarker CSV path                        |
| `--inputfile`      | `configFiles/config_scaffold.txt` | Cell/scaffold config                           |

Shared cluster defaults are in [`scripts/job_defaults.env`](scripts/job_defaults.env). Default submit profile is **gpu** (matches `./scripts/build_drac.sh --cuda`). Use `--profile cpu` for CPU-only builds.

Examples:

```bash
# Default binary at build/bin/testRun
./scripts/submit_testrun.sh your.email@mail.mcgill.ca

# In-tree binary from an older layout
./scripts/submit_testrun.sh your.email@mail.mcgill.ca --testrun bin/testRun

# Preview the sbatch command
./scripts/submit_testrun.sh your.email@mail.mcgill.ca --dry-run
```

Logs go to `logs/testrun_<array>_<jobid>.out` (and `.err`). Simulation artifacts go under **`output/job_<jobid>/`** (biomarker CSV, `tgf_line.csv`, patch dumps, etc.). Each Slurm job also writes **`output/job_<jobid>/run_params.json`** with resolved simulation parameters, the full `testRun` command line, and Slurm metadata. Override output location with `--output-dir` on the submit script or `--output-dir` / `$IVDBM_OUTPUT_DIR` when running `testRun` locally.

## Related docs

- [`configFiles/README.md`](configFiles/README.md) - chemical environment JSON and text config files
- [`src/Diffusion3D/README.md`](src/Diffusion3D/README.md) - Diffusion3D module, tests, and ParaView demo
