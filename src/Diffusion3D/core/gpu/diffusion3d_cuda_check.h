#ifndef DIFFUSION3D_CUDA_CHECK_H
#define DIFFUSION3D_CUDA_CHECK_H

/**
 * @file diffusion3d_cuda_check.h
 * @brief CUDA/cuFFT error checks that abort on failure (safe under `NDEBUG`).
 *
 * Include only from CUDA translation units under `core/gpu/` or from `.cpp`
 * files compiled with `-DDIFFUSION3D_CUDA` (see `DIFFUSION3D_CUDA` CMake option).
 */

#include <cstdio>
#include <cstdlib>

#include <cuda_runtime.h>
#include <cufft.h>

inline void cuda_ok(cudaError_t e, const char *msg) {
  if (e != cudaSuccess) {
    std::fprintf(stderr, "CUDA error: %s (%s)\n", cudaGetErrorString(e), msg);
    std::abort();
  }
}

inline void cufft_ok(cufftResult e, const char *msg) {
  if (e != CUFFT_SUCCESS) {
    std::fprintf(stderr, "cuFFT error: code %d (%s)\n", static_cast<int>(e),
                 msg);
    std::abort();
  }
}

#endif // DIFFUSION3D_CUDA_CHECK_H
