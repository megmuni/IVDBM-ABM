#ifndef GPU_DIFFUSION3D_CUH
#define GPU_DIFFUSION3D_CUH

/**
 * @file gpu_diffusion3d.cuh
 * @brief CUDA translation unit include: `diffusion3d_common.h` plus runtime headers when compiling `.cu`.
 */

#include "diffusion3d_common.h"

#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif

#endif
