/**
 * @file diffusion3d_fft_scratch_stub.cpp
 * @brief Stub FFT scratch API when built without CUDA.
 */

#include <cassert>

#include "diffusion3d_fft_scratch.h"

DiffusionFftScratch *diffusion3d_fft_scratch_create(int, int, int, int, int, int)
{
    return nullptr;
}

void diffusion3d_fft_scratch_destroy(DiffusionFftScratch *) {}

SubstepPlan diffusion3d_fft_scratch_rebuild_operator(DiffusionFftScratch *, double, double, double, double)
{
    assert(false && "diffusion3d_fft_scratch_rebuild_operator requires GPU_DIFFUSE=ON");
    return SubstepPlan{1, 0.0};
}

void diffusion3d_fft_scratch_apply(DiffusionFftScratch *, Diffusion3DContext &)
{
    assert(false && "diffusion3d_fft_scratch_apply requires GPU_DIFFUSE=ON");
}
