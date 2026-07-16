# Diffusion3D

Updates chemical concentration fields on the 3D scaffold grid each tick. The live model uses this for cytokine diffusion (TNF, TGF, IL-1β, O₂). Three backends: CPU stencil, GPU stencil, and GPU FFT.

GPU diffusion uses extra memory for the working grid. The default 3.1 mm world needs an H100 (~5 GB); use 1 mm or smaller on other GPUs.

## Layout

- `core/` - field grids, settings, steppers, VTK export helpers
- `demo/` - standalone ParaView demo (`diffusion3d_demo`)
- Tests: [`../tests/Diffusion3D/`](../tests/Diffusion3D/)

## Build and test

From the repository root:

```bash
cmake -S . -B build -DBUILD_SRC_TESTS=ON
cmake --build build --target diffusion3d_tests
ctest -R diffusion3d --test-dir build
```

GPU tests (CUDA toolkit required):

```bash
cmake -S . -B build \
  -DBUILD_SRC_TESTS=ON \
  -DDIFFUSION3D_CUDA=ON \
  -DDIFFUSION3D_CUDA_ARCHITECTURES=89
cmake --build build --target diffusion3d_tests
ctest -R diffusion3d --test-dir build
```

Use `90` for H100, `89` for Ada (RTX 4050). On DRAC: `./scripts/build_drac.sh --cuda --tests --arch 90`.

Module only (without the full ABM):

```bash
cmake -S src/Diffusion3D -B build -DDIFFUSION3D_TESTS=ON
cmake --build build --target diffusion3d_tests
ctest -R diffusion3d --test-dir build
```

| Flag                             | Default | Purpose                                          |
| -------------------------------- | ------- | ------------------------------------------------ |
| `BUILD_DIFFUSION3D`              | ON      | Build module from root `src/CMakeLists.txt`      |
| `BUILD_SRC_TESTS`                | ON      | Build Catch2 tests under `src/tests/`            |
| `DIFFUSION3D_TESTS`              | ON      | Tests when configuring only `-S src/Diffusion3D` |
| `DIFFUSION3D_DEMO`               | ON      | Build `diffusion3d_demo`                         |
| `DIFFUSION3D_CUDA`               | OFF     | GPU stencil and FFT kernels                      |
| `DIFFUSION3D_CUDA_ARCHITECTURES` | `89`    | GPU arch (`89` Ada, `90` H100, or `native`)      |

## ParaView demo

`diffusion3d_demo` runs all three backends and writes VTK time series:

| Output folder            | Algorithm                       | Needs CUDA |
| ------------------------ | ------------------------------- | ---------- |
| `output/cpu_series/`     | CPU stencil                     | No         |
| `output/gpu_series/`     | GPU stencil                     | Yes        |
| `output/fft_series/`     | GPU FFT                         | Yes        |
| `output/compare_series/` | CPU vs GPU stencil (`abs_diff`) | Yes        |

CPU only:

```bash
cmake -S src/Diffusion3D -B build-demo -DDIFFUSION3D_DEMO=ON
cmake --build build-demo --target diffusion3d_demo
cd build-demo/demo && ./diffusion3d_demo
```

All three backends:

```bash
cmake -S src/Diffusion3D -B build-demo-gpu \
  -DDIFFUSION3D_DEMO=ON -DDIFFUSION3D_CUDA=ON \
  -DDIFFUSION3D_CUDA_ARCHITECTURES=89
cmake --build build-demo-gpu --target diffusion3d_demo
cd build-demo-gpu/demo && ./diffusion3d_demo
```

In ParaView: File --> Open --> `output/cpu_series/diffusion3d_t*.vti` --> Group files as time series.

Demo setup: 8^3 grid, point source at center, 4 macro steps, 5 frames each.

## Export from code

Helpers in [`core/vtk_export.h`](core/vtk_export.h). No VTK library dependency; files are plain XML.

Single species:

```cpp
#include "vtk_export.h"

export_scalar_field_to_vti(grid.grid(species_id), grid.h_,
                           "output/concentration_t5.vti", "concentration");
```

All species in one file:

```cpp
export_multi_species_field_to_vti(multi_grid, "output/all_species_t5.vti");
```

Time series (one file per step):

```cpp
char path[256];
std::snprintf(path, sizeof(path), "output/diffusion3d_t%d.vti", step);
export_multi_species_field_to_vti(grid, path);
```

CPU vs GPU comparison:

```cpp
export_cpu_gpu_compare_to_vti(cpu_u, gpu_u, nx, ny, nz, h,
                              "output/cpu_gpu_compare.vti");
// Arrays: cpu_u, gpu_u, abs_diff
```

VTK spacing follows grid spacing `h`. Use `viz_multiplier` when concentrations are small and ParaView auto-range is hard to read.
