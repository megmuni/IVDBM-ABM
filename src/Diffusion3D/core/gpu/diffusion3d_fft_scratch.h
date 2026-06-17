#ifndef DIFFUSION3D_FFT_SCRATCH_H
#define DIFFUSION3D_FFT_SCRATCH_H

#include <vector>

#include "diffusion3d_timestep.h"

/**
 * @file diffusion3d_fft_scratch.h
 * @brief Opaque GPU workspace for GpuFftPrecomputed (CUDA .cu vs CPU stub .cpp).
 */

struct DiffusionFftScratch;

DiffusionFftScratch *diffusion3d_fft_scratch_create(int nx, int ny, int nz,
                                                    int fft_x, int fft_y, int fft_z);

void diffusion3d_fft_scratch_destroy(DiffusionFftScratch *s);

SubstepPlan diffusion3d_fft_scratch_rebuild_operator(DiffusionFftScratch *s,
                                                     double h,
                                                     double D,
                                                     double tick_dt,
                                                     double safety);

/**
 * @brief Apply the precomputed composed operator to `field` in-place (one convolution).
 *
 * `field.size()` must equal nx*ny*nz stored in the scratch workspace.
 */
void diffusion3d_fft_scratch_apply(DiffusionFftScratch *s, std::vector<double> &field);

#endif
