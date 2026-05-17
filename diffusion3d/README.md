# Standalone 3D Diffusion Demo

Structured-grid 3D heat equation: explicit Euler plus a 6-point face Laplacian (CPU/GPU stencil), or an optional GPU FFT path (`GpuFftPrecomputed`, `GPU_DIFFUSE=ON`) with one padded convolution per macro tick. Demo writes ParaView `.vti`.

## Key files

| File                            | Role                                                           |
| ------------------------------- | -------------------------------------------------------------- |
| `CMakeLists.txt`                | `GPU_DIFFUSE`, Catch2 tests                                    |
| `src/main.cpp`                  | Demo and `.vti` export                                         |
| `src/diffusion3d_common.h`      | `Diffusion3DContext`, `diffusion3d_step_euler_cpu` / `gpu`     |
| `src/diffusion3d_timestep.*`    | `compute_stability_constraint`, `plan_substeps`, `SubstepPlan` |
| `src/diffusion3d_solver.*`      | `DiffusionSolver` (`configure_tick`, `advance_tick`)           |
| `src/diffusion3d_fft_scratch.*` | FFT workspace (cuFFT or stub)                                  |
| `src/gpu_diffusion3d.*`         | CUDA stencil; stub when CUDA off                               |
| `src/convolutionFFT3D.cu`       | FFT helpers (GPU build)                                        |
| `src/vtk_export.h`              | VTK ImageData export                                           |

## Prerequisites

| Dependency          | CPU | GPU |
| ------------------- | --- | --- |
| CMake >= 3.18       | yes | yes |
| C++17               | yes | yes |
| CUDA (nvcc + cuFFT) | —   | yes |

## Main objects

| Object                           | Role                                                                                                                                                                                                                                                                                 |
| -------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `Diffusion3DContext`             | Grid (`nx`, `ny`, `nz`), fields `u` / `u_next`, spacing `h`, diffusivity `D`, GPU block sizes. Index: `(z * ny + y) * nx + x`. Build with `make(nx,ny,nz,h,D)` or `make_with_reference_dt(..., reference_dt)` for an optional `stored_dt` (not used by the solver).                  |
| `DiffusionParams`                | `safety` (see below), `DiffusionBackend`, optional `fft_real_extent_*` (full even R2C real length per axis, or 0 to auto-size).                                                                                                                                                      |
| `DiffusionSolver`                | After `configure_tick(ctx, tick_dt)`, call `advance_tick(ctx)` once per simulation tick to advance `ctx.u` by one macro interval (diffusion only; no sources or reactions). Chooses `n_sub` / `dt_sub` and dispatches stencil or FFT. `resolved_backend()` after `Auto` is resolved. |
| `SubstepPlan`                    | `n_sub`, `dt_sub` from `plan_substeps(tick_dt, dt_max)`.                                                                                                                                                                                                                             |
| `DiffusionFftScratch`            | Opaque GPU workspace; C API `diffusion3d_fft_scratch_*` in `diffusion3d_fft_scratch.h`. Owned by the solver for `GpuFftPrecomputed`.                                                                                                                                                 |
| `DiffusionBackend`               | `CpuStencil`, `GpuStencil`, `GpuFftPrecomputed`, `Auto` (`GpuFftPrecomputed` if `GPU_DIFFUSE=ON`, else `CpuStencil`).                                                                                                                                                                |
| `diffusion3d_step_euler_cpu/gpu` | Single explicit Euler micro-step; pass `dt_sub` explicitly.                                                                                                                                                                                                                          |

FFT vs stencil: padded FFT uses periodic boundaries on the torus; the stencil uses an open box (missing face neighbors omitted). Compare interiors in tests.

## Usage

### Workflow

1. **Grid and physics** — Set `nx`, `ny`, `nz`, `h`, and `D` on the context (via `Diffusion3DContext::make` or `make_with_reference_dt`). These define the discrete Laplacian and the PDE coefficients. Fill `ctx.u` with your initial field.

2. **Macro tick** — Choose `tick_dt`: the real time you want one simulation tick to represent (e.g. one world step of 30 minutes). It may be larger than a single stable explicit step; the solver will split it.

3. **Solver config** — Build `DiffusionParams`: set `backend` (`CpuStencil`, `GpuStencil`, `GpuFftPrecomputed`, or `Auto`), `safety` (see next subsection), and FFT overrides if you use `GpuFftPrecomputed`. Construct `DiffusionSolver solver(params)`.

4. **Configure once per stable setup** — Call `solver.configure_tick(ctx, tick_dt)`. That computes `dt_max = compute_stability_constraint(ctx.h, ctx.D, params.safety)`, plans `(n_sub, dt_sub) = plan_substeps(tick_dt, dt_max)`, and for the FFT backend precomputes the composed operator for this tick. You can log `solver.n_sub()` and `solver.dt_sub()` to see how hard the tick was split.

5. **Advance each tick** — Each simulation step: call `solver.advance_tick(ctx)` once. Do not loop `dt_sub` yourself. Stencil backends run `n_sub` Euler micro-steps; FFT runs one padded convolution per tick.

6. **Reconfigure when inputs change** — If `nx`, `ny`, `nz`, `h`, `D`, `tick_dt`, or `params.safety` / backend / FFT extents change, call `configure_tick` again before the next `advance_tick`. Typical case: `D` or `h` updated from another subsystem (e.g. swelling).

Reference: `src/main.cpp` (export loop, GPU paths behind `GPU_DIFFUSE`). For a host world, map your patch spacing to `h` and your clock step to `tick_dt`.

### Safety (`DiffusionParams::safety`)

`safety` must lie in `(0, 1]`. It multiplies the stability-constraint bound used for the 3D six-point explicit Euler stencil: `compute_stability_constraint(h, D, safety) = safety * h^2 / (6*D)` when `D > 0`. That value is the largest allowed micro-step `dt_sub` before splitting.

- **`safety == 1.0`** (default) — use the bound as written; fewest substeps for a given `tick_dt`.
- **`safety < 1.0`** — shrink `dt_max`, so each micro-step is more conservative: more `n_sub` for the same `tick_dt`, useful if you want margin against boundary effects, coupling to other operators, or distrust of the idealized bound.

The FFT precomputed path uses the same `safety` when composing the kernel so the discrete operator stays consistent with the stencil stability story.

### Units

Keep `h`, `D`, and `tick_dt` in consistent units (e.g. mm, mm^2/min, minutes).

### Symbol quick reference

| Symbol            | Meaning                                                                                                 |
| ----------------- | ------------------------------------------------------------------------------------------------------- |
| `n`               | `nx * ny * nz`                                                                                          |
| `tick_dt`         | Macro interval passed to `configure_tick`                                                               |
| `n_sub`, `dt_sub` | Micro-step count and size; `n_sub * dt_sub == tick_dt` when splitting applies                           |
| `stored_dt`       | Optional on context; not read by `DiffusionSolver`. Use `make` / `make_with_reference_dt` to construct. |

Stability bound for the 3D 6-point stencil: `compute_stability_constraint(h, D, safety)` (when `D > 0`, equals `safety * h^2 / (6*D)`; when `D == 0`, infinite so no diffusion-limited splitting).

## Development

### Build (CPU-only)

```sh
cmake -S . -B build -DGPU_DIFFUSE=OFF
cmake --build build -j
./build/diffusion3d
```

| Folder        | Contents                                               |
| ------------- | ------------------------------------------------------ |
| `cpu_series/` | CPU stencil `.vti`                                     |
| `gpu_series/` | GPU stencil when CUDA on; else repeated initial export |
| `fft_series/` | FFT demo when `GPU_DIFFUSE=ON`                         |

### Build (GPU + FFT)

```sh
cmake -S . -B build-gpu -DGPU_DIFFUSE=ON
cmake --build build-gpu -j
./build-gpu/diffusion3d
```

Set `CUDACXX` or ensure `nvcc` is on `PATH` if needed.

### Visualization (ParaView)

Open `cpu_series/diffusion3d_t*.vti` (and `gpu_series/`, `fft_series/` as built). Index `k` matches the demo loop in `main.cpp`.

## Tests (CTest)

Catch2 via `FetchContent` when `DIFFUSION3D_TESTS=ON` (default).

```sh
cmake -S . -B build -DGPU_DIFFUSE=OFF -DDIFFUSION3D_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

`-DDIFFUSION3D_TESTS=OFF` disables tests.

## Logging

```sh
cmake -S . -B build -DLOGGING=ON
```

Enables `STEP_LOG` and related tracing.
