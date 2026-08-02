# Configuration files

## Simulation config (required)

Scaffold seeding, chemistry, and biology parameters are read from a single JSON file with `world_init`, `chemistry`, and `biology` sections. The C++ loaders extract each section independently.

### Setup

1. Copy the template:
   ```bash
   cp configFiles/simulation_config.template.json configFiles/simulation_config.json
   ```
2. Edit `simulation_config.json`.
3. Run from the project root (default path), or pass a custom file:
   ```bash
   ./build/bin/testRun --config /path/to/simulation_config.json
   ```

### Root fields

| Section      | Purpose                                                     |
| ------------ | ------------------------------------------------------------ |
| `world_init` | Scaffold seeding and alginate composition (see below)       |
| `chemistry`  | Diffusion, baselines, tick length (see below)                |
| `biology`    | Cell rules and hydrogel calibration (see below)              |

### `chemistry` section

| Field                           | Purpose                                                                                                                                         |
| ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `tick_interval_minutes`         | Length of one tick (minutes); template default is 30                                                                                            |
| `channels.count`                | Number of float grids allocated                                                                                                                 |
| `channels.chemotaxis`           | Channel index used for cell chemotaxis (`pcellgrad` = 6 in IVDBM)                                                                               |
| `channels.channel_names`        | Optional labels for documentation                                                                                                               |
| `merge.chemotaxis_from_species` | Species name copied into the chemotaxis field each merge (default `"TGF"`)                                                                      |
| `baseline_total_mass`           | Initial total mass per species; divided over CaAlg patches at startup                                                                           |
| `species[]`                     | Each diffusing species: `id`, `name`, `base_diffusivity_mm2_per_min`, `concentration_channel`, `diffused_channel`, optional `diffusivity_model` |

Effective diffusivity at runtime is per species (optional `diffusivity_model`; default `swelling_ratio`):

| Model                      | Formula                                                                        |
| -------------------------- | ------------------------------------------------------------------------------ |
| `swelling_ratio` (default) | `base_diffusivity_mm2_per_min × Q`                                             |
| `logarithmic_stiffness`    | `slope × ln(E) + intercept` (E in kPa; O₂ uses slope −0.002, intercept 0.0218) |

IVDBM template defaults:

| Species | Base D (mm²/min) | p ch | d ch |
| ------- | ---------------- | ---- | ---- |
| TNF     | 0.00018          | 0    | 4    |
| TGF     | 0.000156         | 1    | 5    |
| IL-1β   | 0.00018          | 2    | 6    |
| O₂      | 0.02172          | 3    | 7    |

### `biology` section

| Group         | Purpose                                                                       |
| ------------- | ----------------------------------------------------------------------------- |
| `cell`        | Shared proliferation timing and cytokine synthesis feedback (k2–k17)          |
| `stem`        | MSC OCR, apoptosis, migration, cytokines, ECM, proliferation, differentiation |
| `progen`      | Pre-NP OCR, apoptosis, migration, cytokines, aggrecan                         |
| `np`          | NP OCR, apoptosis, migration, collagen/aggrecan synthesis vs time             |
| `biomaterial` | Ca-Alg elastic modulus, pore size, mass loss, swelling ratio (c1–c18)         |

Loaded biology values are printed to stdout at startup. When `IVDBM_RUN_PARAMS_JSON` is set, they are also written into `run_params.json`.

### `world_init` section

Scaffold seeding and alginate hydrogel composition:

| Field                        | Purpose                                                          |
| ----------------------------- | ----------------------------------------------------------------- |
| `msc_count`                  | Initial MSC seed count; `0` seeds at the default density (1e6 cells/mL) |
| `alginate.wv_percent`        | Alginate concentration (% w/v)                                   |
| `alginate.high_mw_ratio`     | Ratio component of high-MW alginate                              |
| `alginate.low_mw_ratio`      | Ratio component of low-MW alginate                                |
| `alginate.ca_mm`             | Calcium crosslinker concentration (mM)                            |
| `peptide`                    | Peptide conjugation type; only used when compiled with `PEPTIDE_BM` |

---

### Which file controls what

| Concern                                 | File / section                         |
| --------------------------------------- | -------------------------------------- |
| Diffusion, baselines, tick length       | `simulation_config.json` → `chemistry` |
| Initial cell counts, alginate/Ca inputs | `simulation_config.json` → `world_init` |
| Cell rules and hydrogel calibration     | `simulation_config.json` → `biology`   |
