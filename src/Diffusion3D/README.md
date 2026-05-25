# Diffusion3D Module — Phase I

This module implements multi-species scalar field diffusion with unified double-precision storage and shared pointer support. It is designed for integration with ABM frameworks and supports extensible PDE stepping strategies.

## Structure

- `core/` - Core diffusion engine, field grids, settings, stepper implementations, and VTK export
- `demo/` - ParaView `.vti` time-series demo (`diffusion3d_demo`) for all three algorithms
- `tests/` - Unit and integration tests (support helpers under `tests/support/`)

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
| `DIFFUSION3D_DEMO` | ON | Build ParaView demo (`diffusion3d_demo`) |
| `DIFFUSION3D_CUDA` | OFF | Enable GPU/FFT kernels in `core/gpu/` and GPU steppers |
| `DIFFUSION3D_CUDA_ARCHITECTURES` | `89` | CUDA arch for GPU targets (`89` = sm_89, Ada/RTX 4050; or `native`, `89;90`) |

Tests use `MultiSpeciesFieldGrid`, `MultiSpeciesDiffusionEngine`, and `MakeMultiSpeciesDiffusionEngine()` exclusively — no legacy solver or context types remain in the module.

CUDA build example (RTX 4050 / sm_89):

```bash
cmake -S src/Diffusion3D -B build-cuda \
  -DDIFFUSION3D_TESTS=ON \
  -DDIFFUSION3D_CUDA=ON \
  -DDIFFUSION3D_CUDA_ARCHITECTURES=89
cmake --build build-cuda --target diffusion3d_tests
ctest -R diffusion3d --test-dir build-cuda
```

## ParaView demo (`diffusion3d_demo`)

Runs all three diffusion backends and writes VTK ImageData time series for ParaView:

| Output folder | Algorithm | Requires CUDA |
| ------------- | --------- | ------------- |
| `output/cpu_series/` | CPU explicit Euler stencil | no |
| `output/gpu_series/` | GPU 6-point stencil | yes |
| `output/fft_series/` | GPU FFT precomputed | yes |
| `output/compare_series/` | CPU vs GPU stencil at final step (`abs_diff` array) | yes |

Build and run from the module directory:

```bash
# CPU only (cpu_series/)
cmake -S src/Diffusion3D -B build-demo -DDIFFUSION3D_DEMO=ON
cmake --build build-demo --target diffusion3d_demo
cd build-demo/demo && ./diffusion3d_demo

# All three algorithms
cmake -S src/Diffusion3D -B build-demo-gpu \
  -DDIFFUSION3D_DEMO=ON -DDIFFUSION3D_CUDA=ON \
  -DDIFFUSION3D_CUDA_ARCHITECTURES=89
cmake --build build-demo-gpu --target diffusion3d_demo
cd build-demo-gpu/demo && ./diffusion3d_demo
```

In ParaView: **File → Open** → select `output/cpu_series/diffusion3d_t*.vti` → enable **Group files as time series**.

Scenario matches the legacy standalone demo: 8³ grid, point source at center, `tick_dt=2`, `D=0.1`, `h=1`, 4 macro steps (5 frames each).

## ParaView visualization (API)

Export simulation fields as VTK **ImageData** (`.vti`) and open them in [ParaView](https://www.paraview.org/). Header-only helpers live in [`core/vtk_export.h`](core/vtk_export.h).

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
| **Demo** | Run `diffusion3d_demo` — see [ParaView demo](#paraview-demo-diffusion3d_demo) above |

No VTK library dependency — files are plain XML written with the standard library.
