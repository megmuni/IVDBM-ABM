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
cp configFiles/simulation_config.template.json configFiles/simulation_config.json
mkdir -p output
```

Chemistry field details are in [configFiles/README.md](configFiles/README.md).

## What to change

| What you want to change                         | Where to set it                                     | When                  |
| ----------------------------------------------- | --------------------------------------------------- | --------------------- |
| Cytokine diffusivities, baselines, tick length  | `configFiles/simulation_config.json` → `chemistry`  | Before a run          |
| Cell seeding, alginate, scaffold layout         | `configFiles/simulation_config.json` → `world_init` | Before a run          |
| Cell rule and hydrogel calibration coefficients | `configFiles/simulation_config.json` → `biology`    | Before a run          |
| Tick count, world size, paths, output           | `testRun` flags or `submit_testrun.sh` flags        | Each run              |
| Shared DRAC simulation defaults                 | `scripts/job_defaults.env`                          | Optional              |
| Model type, 3D, GPU diffusion, biomarker output | `src/common.h`                                      | Rebuild after editing |
| Job time limit, CPUs, GPU choice                | `submit_testrun.sh` flags                           | Each DRAC job         |
| Unit test job defaults                          | `scripts/test_defaults.env` or `submit_tests.sh`    | Test jobs only        |

Command-line flags override `job_defaults.env`, which overrides the built-in defaults. With no flags on your own machine: 432 ticks, 3.1 mm cube. On DRAC with no extra flags: 24 ticks, 1 mm cube.

Biology settings are read from JSON when the run starts. Loaded biology parameter values are printed to stdout and recorded in `run_params.json` when `IVDBM_RUN_PARAMS_JSON` is set.

### Changing a biology/hydrogel parameter value

Most calibration work is just editing numbers. Open `configFiles/simulation_config.json`, find the field under `biology` (see [configFiles/README.md](configFiles/README.md) for the full field list), and change its value:

```jsonc
"stem": {
  "proliferation": {
    "tgf_threshold": 10   // <-- edit this
  }
}
```

No rebuild needed, just rerun `testRun`.

### Adding a brand-new biology parameter

If you need a parameter that doesn't exist yet (not just a new value for an existing one), it lives in two places that must stay in sync: a C++ struct that mirrors the JSON, and a small static array the cell/hydrogel formulas actually read from. Using "add a new stem cell proliferation coefficient" as an example:

1. Add the field to the matching struct in `src/Agent/biology_parameters_config.h`:

   ```cpp
   struct StemProliferationParams {
     double tgf_threshold = 0;
     double tnf_effect = 0;
     double il1beta_effect = 0;
     double elasticity_effect = 0;
     double my_new_effect = 0;   // <-- new field
   };
   NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StemProliferationParams, tgf_threshold,
                                      tnf_effect, il1beta_effect,
                                      elasticity_effect, my_new_effect)  // <-- add it here too
   ```

2. Add a matching named slot in `src/Agent/Usr_Agents/Cell.h` (the array resizes itself, and a compile-time check will fail loudly if this step and step 1 ever get out of sync):

   ```cpp
   enum ProliferationIdx {
     PROLIFERATION_TGF_THRESHOLD = 0,
     PROLIFERATION_TNF_EFFECT,
     PROLIFERATION_IL1BETA_EFFECT,
     PROLIFERATION_ELASTICITY_EFFECT,
     PROLIFERATION_MY_NEW_EFFECT,   // <-- new slot
     PROLIFERATION_COUNT
   };
   ```

3. Wire the JSON value into that slot in `apply_biology_parameters()` in `src/Agent/biology_parameters_config.cpp`:

   ```cpp
   Stem::proliferation[Stem::PROLIFERATION_MY_NEW_EFFECT] =
       static_cast<float>(st.proliferation.my_new_effect);
   ```

4. Add the key to `configFiles/simulation_config.template.json` (and your local `simulation_config.json`), under `biology.stem.proliferation.my_new_effect`.

5. Use it in a formula, e.g. in `src/Agent/Usr_Agents/Cell.cpp`, via `Stem::proliferation[Stem::PROLIFERATION_MY_NEW_EFFECT]`.

Rebuild after these changes. Reordering or renaming an existing slot works the same way - just edit the enum name/position in step 2 and its one reference in step 3; formulas elsewhere keep working since they refer to the name, not a position.

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

| Parameter         | Flag or file              | Notes                                                                                                         |
| ----------------- | ------------------------- | ------------------------------------------------------------------------------------------------------------- |
| Tick count        | `--numticks`              | 24 ≈ 12 h, 48 ≈ 1 day, 432 ≈ 9 days                                                                           |
| World size (mm)   | `--wxw`, `--wyw`, `--wzw` | Local default 3.1 mm; DRAC default 1 mm. GPU: 3.1 mm needs an H100 (~5 GB); use 1 mm or smaller on other GPUs |
| Simulation config | `--config`                | `configFiles/simulation_config.json` — world_init + chemistry + biology sections                              |
| Output directory  | `--output-dir`            | All CSVs; on DRAC also `run_params.json`                                                                      |
| Biomarker CSV     | `--outputfile`            | Default `<output-dir>/Output_Biomarkers.csv`                                                                  |
| ParaView export   | `--paraview`              | Writes `<output-dir>/paraview/patches_t*.vti` and `chem_t*.vti` every 12 ticks                                |

Local:

```bash
# --- simulation ---
NUMTICKS=48
WXW=3.1
WYW=3.1
WZW=3.1
SIMULATION_CONFIG="configFiles/simulation_config.json"
OUTPUT_DIR="output/local_${NUMTICKS}ticks"
OUTPUTFILE="${OUTPUT_DIR}/Output_Biomarkers.csv"

./build/bin/testRun \
  --output-dir "${OUTPUT_DIR}" \
  --numticks "${NUMTICKS}" \
  --wxw "${WXW}" --wyw "${WYW}" --wzw "${WZW}" \
  --config "${SIMULATION_CONFIG}" \
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
SIMULATION_CONFIG="configFiles/simulation_config.json"
OUTPUTFILE=""   # leave empty for default: <output-dir>/Output_Biomarkers.csv

# --- cluster ---
PROFILE="gpu"              # gpu after ./scripts/build_drac.sh --cuda; cpu for CPU-only build
WALLTIME="0-02:00:00"      # raise for long runs (submit default is 5 minutes)
CPUS=32
MEM="32000M"

SUBMIT=(./scripts/submit_testrun.sh "${EMAIL}" \
  --numticks "${NUMTICKS}" \
  --wxw "${WXW}" --wyw "${WYW}" --wzw "${WZW}" \
  --config "${SIMULATION_CONFIG}" \
  --profile "${PROFILE}" \
  --time "${WALLTIME}" \
  --cpus "${CPUS}" \
  --mem "${MEM}")

[[ -n "${OUTPUTFILE}" ]] && SUBMIT+=(--outputfile "${OUTPUTFILE}")

"${SUBMIT[@]}"
```

Append `--dry-run` to preview without submitting.

`./scripts/submit_testrun.sh "${EMAIL}" 12` is the same as `--numticks 12`.

To change shared DRAC simulation defaults, edit `scripts/job_defaults.env`. To change biology, edit files under `configFiles/` (see [configFiles/README.md](configFiles/README.md)).

## ParaView visualization

Pass `--paraview` to `testRun` to export VTK ImageData (`.vti`) time series under `<output-dir>/paraview/`:

| File                  | Contents                                       |
| --------------------- | ---------------------------------------------- |
| `patches_t<tick>.vti` | Patch occupancy (`occupied`: 0/1)              |
| `chem_t<tick>.vti`    | Species concentrations (TNF, TGF, IL1beta, o2) |

Exports occur every 12 ticks (tick 0, 12, 24, …). In ParaView: **File --> Open** --> select `chem_t*.vti` or `patches_t*.vti` --> enable **Group files as time series**.

Example:

```bash
./build/bin/testRun \
  --numticks 24 \
  --paraview \
  --wxw 1 --wyw 1 --wzw 1 \
  --config configFiles/simulation_config.template.json \
  --output-dir output/paraview_smoke
```

Diffusion3D-only ParaView demos are documented in [src/Diffusion3D/README.md](src/Diffusion3D/README.md).

## Results

| What                      | Local run                               | DRAC job                                      |
| ------------------------- | --------------------------------------- | --------------------------------------------- |
| Progress and errors       | Terminal                                | `logs/testrun_<array>_<jobid>.out` and `.err` |
| CSVs and calibration data | `<output-dir>/` (default `output/`)     | `output/job_<jobid>/`                         |
| Run settings snapshot     | Only if `$IVDBM_RUN_PARAMS_JSON` is set | `<output-dir>/run_params.json`                |

Main output files (default scaffold build):

| File                                       | When       | Contents                                        |
| ------------------------------------------ | ---------- | ----------------------------------------------- |
| `Output_Biomarkers.csv`                    | Each tick  | Cell counts, cytokines, ECM, scaffold mechanics |
| `tgf_line.csv`                             | Each tick  | TGF concentration along the x-axis              |
| `SensitivityAnalysis/FinalTotalChemVR.dat` | End of run | Chemical totals for calibration                 |

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

| `--suite`     | What it runs           | GPU needed                                          |
| ------------- | ---------------------- | --------------------------------------------------- |
| `chemistry`   | Chemistry tests        | No                                                  |
| `diffusion3d` | Diffusion3D tests      | Mixed                                               |
| `world`       | World tests            | No                                                  |
| `cpu`         | All CPU-labelled tests | No                                                  |
| `gpu`         | GPU-labelled tests     | Yes                                                 |
| `all`         | Full suite             | Use CUDA build and `--profile gpu` for GPU coverage |

Test job logs: `logs/tests_<array>_<jobid>.out` / `.err`. CTest output: `output/tests_<jobid>/ctest.log`.

Shared test job defaults: `scripts/test_defaults.env`.

## Related docs

- [configFiles/README.md](configFiles/README.md): chemistry JSON and cell/scaffold text configs
- [src/Diffusion3D/README.md](src/Diffusion3D/README.md): diffusion module, tests, ParaView demo
