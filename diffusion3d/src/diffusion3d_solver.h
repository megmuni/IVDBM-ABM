#ifndef DIFFUSION3D_SOLVER_H
#define DIFFUSION3D_SOLVER_H

#include <cassert>
#include <cstddef>

#include "diffusion3d_common.h"

/**
 * @file diffusion3d_solver.h
 * @brief High-level diffusion advance: substep planning plus CPU, GPU stencil, or GPU FFT backend (CUDA-free header).
 */

struct DiffusionFftScratch;

/**
 * @brief Which hardware path advances the field for one macro interval `tick_dt`.
 *
 * Auto picks GpuFftPrecomputed when built with GPU_DIFFUSE, otherwise CpuStencil.
 * GpuFftPrecomputed is CUDA-only: one padded FFT per advance with spectrum (K_hat)^N (N = substep count).
 */
enum class DiffusionBackend
{
    CpuStencil,
    GpuStencil,
    /** @brief Padded cuFFT conv with combined kernel spectrum; periodic BC vs open stencil. */
    GpuFftPrecomputed,
    Auto,
};

/**
 * @brief User-tunable solver configuration (does not duplicate grid scalars in Diffusion3DContext).
 *
 * safety scales the explicit-Euler stability bound (see diffusion3d_timestep).
 */
struct DiffusionParams
{
    double safety = 1.0;
    DiffusionBackend backend = DiffusionBackend::Auto;
    /**
     * Full padded real R2C grid length per axis (even). 0 = derive automatically.
     * Not a delta in voxels; see README (Usage / main objects).
     */
    int fft_real_extent_x = 0;
    int fft_real_extent_y = 0;
    int fft_real_extent_z = 0;
};

/**
 * @brief Owns optional GPU device buffers; advances `ctx.u` by `tick_dt` via stable substeps.
 *
 * GpuFftPrecomputed stores an opaque `DiffusionFftScratch` (cuFFT + padded buffers).
 */
class DiffusionSolver
{
public:
    explicit DiffusionSolver(DiffusionParams params);

    DiffusionSolver(const DiffusionSolver &) = delete;
    DiffusionSolver &operator=(const DiffusionSolver &) = delete;

    ~DiffusionSolver();

    /** @brief Concrete backend: `params.backend` with `Auto` expanded (never `Auto` here). */
    DiffusionBackend resolved_backend() const { return resolved_backend_; }

    /**
     * @brief Configure a fixed macro tick. Must be called before `advance_tick`.
     *
     * This computes the CFL substep plan (n_sub, dt_sub) for the current `ctx` parameters and stores it.
     * FFT backend also precomputes the composed kernel spectrum for this tick.
     */
    void configure_tick(Diffusion3DContext &ctx, double tick_dt);

    /**
     * @brief Advance exactly one configured macro tick (sole high-level per-tick API after `configure_tick`).
     *
     * CpuStencil/GpuStencil: runs `n_sub` Euler micro-steps at `dt_sub`.
     * GpuFftPrecomputed: one pad/R2C/modulate/C2R/unpad using precomputed spectrum (K_hat_one_step)^n_sub.
     */
    void advance_tick(Diffusion3DContext &ctx);

    /** @brief Configured macro tick. Requires `configure_tick` to have been called. */
    double tick_dt() const { return tick_dt_; }
    int n_sub() const { return n_sub_; }
    double dt_sub() const { return dt_sub_; }

private:
    void ensure_gpu_buffers(std::size_t n);
    void free_gpu_buffers();
    void ensure_fft_scratch(Diffusion3DContext &ctx);

    DiffusionParams params_;
    /** @brief Same as `params_.backend` except `Auto` is replaced at construction. */
    DiffusionBackend resolved_backend_;

    bool configured_{false};
    double tick_dt_{0.0};
    int n_sub_{1};
    double dt_sub_{0.0};
    int cfg_nx_{0};
    int cfg_ny_{0};
    int cfg_nz_{0};
    double cfg_h_{0.0};
    double cfg_D_{0.0};

    double *d_u_{nullptr};
    double *d_u_next_{nullptr};
    std::size_t gpu_n_{0};

    DiffusionFftScratch *fft_scratch_{nullptr};
    int fft_cached_nx_{-1};
    int fft_cached_ny_{-1};
    int fft_cached_nz_{-1};
    int fft_cached_fx_{-1};
    int fft_cached_fy_{-1};
    int fft_cached_fz_{-1};
};

#endif // DIFFUSION3D_SOLVER_H
