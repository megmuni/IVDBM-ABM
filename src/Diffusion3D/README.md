# Diffusion3D Module — Phase I

This module implements multi-species scalar field diffusion with unified double-precision storage and shared pointer support. It is designed for integration with ABM frameworks and supports extensible PDE stepping strategies.

## Structure

- `core/` - Core diffusion engine, field grids, settings, stepper implementations, and VTK export
- `tests/` - Unit and integration tests

## Key Types

- `MultiSpeciesFieldGrid` - Container for per-species 3D scalar fields on shared domain
- `MultiSpeciesDiffusionSettings` - Per-species diffusivity, CFL safety, and algorithm selection
- `MultiSpeciesDiffusionEngine` - Abstract strategy interface for PDE stepping
- `ExplicitMultiSpeciesHeatStepper` - CPU explicit Euler (Phase I reference implementation)
- `MakeMultiSpeciesDiffusionEngine()` - Factory for engine creation with validation
- `DiffusionAlgorithm` - Enum: ExplicitHeatEquation, GpuStencil, GpuFftPrecomputed
- `vtk_export.h` - ParaView `.vti` export for scalar and multi-species fields

## Build

From the repo root:

```bash
cmake -S . -B build -DBUILD_DIFFUSION3D=ON
cmake --build build --target diffusion3d_tests
ctest -R diffusion3d --test-dir build
```

Options:

| Flag | Default | Purpose |
| ---- | ------- | ------- |
| `BUILD_DIFFUSION3D` | ON | Build module and tests from root `src/CMakeLists.txt` |
| `DIFFUSION3D_TESTS` | ON | Build Catch2 test executable |
| `DIFFUSION3D_CUDA` | OFF | Enable GPU paths in legacy numerical tests |

Tests live under `tests/` and use Catch2 v2 via `FetchContent`. Legacy kernel sources remain in `diffusion3d/src/` until Phase I migration completes.

## ParaView visualization

Export simulation fields as VTK **ImageData** (`.vti`) and open them in [ParaView](https://www.paraview.org/). Header-only helpers live in [`core/vtk_export.h`](core/vtk_export.h) (adapted from the standalone `diffusion3d` demo).

### Single-species field

```cpp
#include "vtk_export.h"

// After advance_species_interval(...)
export_scalar_field_to_vti(grid.grid(species_id), grid.h_,
                           "output/concentration_t5.vti", "concentration");
```

### Multi-species (one file, one array per species)

```cpp
export_multi_species_field_to_vti(multi_grid, "output/all_species_t5.vti");
// ParaView: color by species_0, species_1, ...
```

### Time series

Write one `.vti` per macro tick (same pattern as the standalone demo):

```cpp
char path[256];
std::snprintf(path, sizeof(path), "output/diffusion3d_t%d.vti", step);
export_multi_species_field_to_vti(grid, path);
```

Then load the sequence in ParaView: **File → Open** → select all `diffusion3d_t*.vti` → **OK** → enable **Group files as time series**.

### CPU vs GPU comparison

```cpp
export_cpu_gpu_compare_to_vti(cpu_u, gpu_u, nx, ny, nz, h,
                              "output/cpu_gpu_compare.vti");
// Arrays: cpu_u, gpu_u, abs_diff
```

### Notes

| Topic | Detail |
| ----- | ------ |
| **Spacing** | VTK `Spacing` is set from grid spacing `h` (physical units consistent with the simulation) |
| **Layout** | `ScalarFieldGrid::data_` uses the same point order as VTK ImageData (`idx = x + nx*(y + ny*z)`) |
| **Scaling** | Use `viz_multiplier` on `export_scalar_field_to_vti` / `export_field_to_vti` when concentrations are small and ParaView auto-range is hard to read |
| **Standalone demo** | The reference time-series workflow is in [`diffusion3d/src/main.cpp`](../../diffusion3d/src/main.cpp); it will link `diffusion3d_core` after I.8 migration |

No VTK library dependency — files are plain XML written with the standard library.
