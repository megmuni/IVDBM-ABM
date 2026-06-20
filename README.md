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
  --numticks 288 \
  --inputfile configFiles/config_scaffold.txt \
  --outputfile output/Output_Biomarkers.csv \
  --wxw 1 --wyw 1 --wzw 1
```

Useful CLI flags (see `./build/bin/testRun --help`):

| Flag                      | Description                                                                            |
| ------------------------- | -------------------------------------------------------------------------------------- |
| `--numticks`              | Number of simulation ticks (1 tick = 30 min by default)                                |
| `--wxw`, `--wyw`, `--wzw` | World width, length, height (mm)                                                       |
| `--inputfile`             | Text config for cells, alginate, and scaffold parameters                               |
| `--outputfile`            | CSV biomarker output path                                                              |
| `--chem-config`           | Path to `chemical_environment.json` (default: `configFiles/chemical_environment.json`) |

With `MODEL_SCAFFOLD` defined in `src/common.h`, defaults match a small alginate scaffold (`configFiles/config_scaffold.txt`, 3.1 mm cube unless overridden on the command line).

## Run `testRun` on DRAC (Slurm)

Use [`scripts/submit_testrun.sh`](scripts/submit_testrun.sh) from the repository root. It calls [`scripts/testrun.sbatch`](scripts/testrun.sbatch) with the correct `#SBATCH` options, your email, and paths to the binary and config files.

```bash
./scripts/submit_testrun.sh your.email@mail.mcgill.ca
```

Common options:

| Flag             | Default                           | Purpose                                        |
| ---------------- | --------------------------------- | ---------------------------------------------- |
| `--testrun PATH` | `build/bin/testRun`               | Executable (relative to repo root or absolute) |
| `--array`        | `0-4`                             | Job array range                                |
| `--account`      | `def-nicoleli`                    | Slurm allocation                               |
| `--numticks`     | `288`                             | Simulation length                              |
| `--inputfile`    | `configFiles/config_scaffold.txt` | Cell/scaffold config                           |

Examples:

```bash
# Default binary at build/bin/testRun
./scripts/submit_testrun.sh your.email@mail.mcgill.ca

# In-tree binary from an older layout
./scripts/submit_testrun.sh your.email@mail.mcgill.ca --testrun bin/testRun

# Preview the sbatch command
./scripts/submit_testrun.sh your.email@mail.mcgill.ca --dry-run
```

Logs go to `logs/testrun_<array>_<jobid>.out`; biomarker CSVs to `output/Output_Biomarkers_<task>.csv`.

## Related docs

- [`configFiles/README.md`](configFiles/README.md) - chemical environment JSON and text config files
- [`src/Diffusion3D/README.md`](src/Diffusion3D/README.md) - Diffusion3D module, tests, and ParaView demo
