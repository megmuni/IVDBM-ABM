# Configuration files

## Chemical environment (required)

Biological values for the ABM chemistry system are read **only** from JSON. C++ does not embed IVDBM diffusivities, baselines, or tick length.

### Setup

1. Copy the template:
   ```bash
   cp configFiles/chemical_environment.template.json configFiles/chemical_environment.json
   ```
2. Edit `chemical_environment.json`.
3. Run from the **project root** (default path), or pass a custom file:
   ```bash
   ./bin/Nyvonna --chem-config /path/to/chemical_environment.json
   ```

### JSON fields (schema version 1)

| Section                         | Purpose                                                                                                           |
| ------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `schema_version`                | Must be `1`                                                                                                       |
| `tick_interval_minutes`         | Length of one ABM tick (minutes); template default is **30** (matches ABM clock)                                  |
| `channels.count`                | Number of float grids allocated                                                                                   |
| `channels.chemotaxis`           | Channel index used for cell chemotaxis (`pcellgrad` = 6 in IVDBM)                                                 |
| `channels.channel_names`        | Optional labels for documentation                                                                                 |
| `merge.chemotaxis_from_species` | Species name copied into the chemotaxis field each merge (default `"TGF"`)                                        |
| `baseline_total_mass`           | Initial total mass per species; divided over CaAlg patches at startup                                             |
| `species[]`                     | Each diffusing species: `id`, `name`, `base_diffusivity_mm2_per_min`, `concentration_channel`, `diffused_channel`, optional `diffusivity_model` |

Effective diffusivity at runtime is **per species** (optional `diffusivity_model`; default `swelling_ratio`):

| Model | Formula |
| ----- | ------- |
| `swelling_ratio` (default) | `base_diffusivity_mm2_per_min × Q` |
| `logarithmic_stiffness` | `slope × ln(E) + intercept` (E in kPa; O₂ uses slope −0.002, intercept 0.0218) |

### IVDBM template defaults

| Species | Base D (mm²/min) | p ch | d ch |
| ------- | ---------------- | ---- | ---- |
| TNF     | 0.00018          | 0    | 3    |
| TGF     | 0.000156         | 1    | 4    |
| IL-1β   | 0.00018          | 2    | 5    |
| O₂      | 0.02172          | 3    | 7    |

---

## World / cells (`config.txt` and variants)

`config.txt`, `config_Scaffold_GH10.txt`, and related files still supply:

- Initial cell counts
- Alginate / Ca crosslinker inputs (scaffold builds)

**Not used for chemistry after JSON load:** `Baseline_TNF`, `Baseline_TGF`, `Baseline_IL1beta`, and legacy `D_*` / `HL` tags in the text file. Set those in `chemical_environment.json` instead. TODO: Need to be removed.
