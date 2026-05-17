#include <cassert>

#include "diffusion3d_common.h"

/**
 * @file gpu_diffusion3d_stub.cpp
 * @brief Stub `diffusion3d_step_euler_gpu` when the project is built without CUDA.
 */

void diffusion3d_step_euler_gpu(Diffusion3DContext &, double, double *&, double *&)
{
    assert(false && "diffusion3d_step_euler_gpu requires GPU_DIFFUSE build");
}