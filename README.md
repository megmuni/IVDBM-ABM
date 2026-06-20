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

To enable GPU diffusion kernels (required on CUDA Slurm nodes for the GPU/FFT chemistry path):

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

## Run `testRun` as a Slurm job

[`slurm_job.sh`](slurm_job.sh) submits a job array on a Compute Canada–style cluster (`def-nicoleli` account, H100 GPUs). Edit the `#SBATCH` lines (account, time, memory, mail user, array range) before submitting.

The script expects the executable at `./bin/testRun` relative to the repo root. After an out-of-tree build, either symlink or copy the binary:

```bash
ln -sf build/bin/testRun bin/testRun
```

Then submit from the repository root:

```bash
module load cuda    # loaded inside slurm_job.sh as well
sbatch slurm_job.sh
```

The default script body runs five array tasks (`--array=0-4`), writing Slurm logs to `test_<task_id>.out` and biomarker CSVs to `output/Output_Biomarkers_<task_id>.csv`. Uncomment the plain `./bin/testRun ...` line (and comment out the GDB line) for production runs without the debugger.

Example production command (already present, commented, in `slurm_job.sh`):

```bash
./bin/testRun \
  --numticks 200 \
  --inputfile config_scaffold.txt \
  --wxw 3 --wyw 3 --wzw 3
```

On the compute node, load CUDA before building or running if GPU diffusion is enabled (`module load cuda`), and match `-DDIFFUSION3D_CUDA_ARCHITECTURES` to the node GPU (e.g. `90` for H100).

## Related docs

- [`configFiles/README.md`](configFiles/README.md) - chemical environment JSON and text config files
- [`src/Diffusion3D/README.md`](src/Diffusion3D/README.md) - Diffusion3D module, tests, and ParaView demo
