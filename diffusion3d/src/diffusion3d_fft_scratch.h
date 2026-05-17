#ifndef DIFFUSION3D_FFT_SCRATCH_H
#define DIFFUSION3D_FFT_SCRATCH_H

#include "diffusion3d_common.h"
#include "diffusion3d_timestep.h"

/**
 * @file diffusion3d_fft_scratch.h
 * @brief Opaque GPU workspace for GpuFftPrecomputed (CUDA .cu vs CPU stub .cpp).
 */

struct DiffusionFftScratch;

/** @brief Allocates device buffers and cuFFT plans (fft_* even, fft_* >= domain + 3x3x3 kernel halo). */
DiffusionFftScratch *diffusion3d_fft_scratch_create(int nx, int ny, int nz,
                                                    int fft_x, int fft_y, int fft_z);

void diffusion3d_fft_scratch_destroy(DiffusionFftScratch *s);

/**
 * @brief Precompute composed operator for a fixed macro tick.
 *
 * Computes `plan_substeps(tick_dt, compute_stability_constraint(h,D,safety))`, then stores the composed
 * kernel spectrum (K_hat_one_step)^n_sub inside `s`.
 *
 * Returns the (n_sub, dt_sub) plan used for the composition.
 *
 * Periodic wrap from the padded torus differs from the open-box stencil at boundaries; compare interiors in tests.
 */
SubstepPlan diffusion3d_fft_scratch_rebuild_operator(DiffusionFftScratch *s,
                                                     double h,
                                                     double D,
                                                     double tick_dt,
                                                     double safety);

/**
 * @brief Apply the precomputed composed operator to `ctx.u` (one convolution).
 *
 * Requires `diffusion3d_fft_scratch_rebuild_operator` to have been called.
 */
void diffusion3d_fft_scratch_apply(DiffusionFftScratch *s, Diffusion3DContext &ctx);

#endif
