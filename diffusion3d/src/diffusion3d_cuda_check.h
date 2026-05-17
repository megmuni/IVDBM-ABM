#ifndef DIFFUSION3D_CUDA_CHECK_H
#define DIFFUSION3D_CUDA_CHECK_H

/**
 * @file diffusion3d_cuda_check.h
 * @brief Always-executed CUDA/cuFFT status checks; safe under `NDEBUG`.
 *
 * Only available in GPU builds/
 */

#include <cstdio>
#include <cstdlib>

#ifdef GPU_DIFFUSE
#include <cuda_runtime.h>
#include <cufft.h>

#endif

#ifdef GPU_DIFFUSE
/** @brief Abort with a CUDA error string if `e` is not `cudaSuccess`. */
inline void cuda_ok(cudaError_t e, const char *msg) {
  if (e != cudaSuccess) {
    std::fprintf(stderr, "CUDA error: %s (%s)\n", cudaGetErrorString(e), msg);
    std::abort();
  }
}

/** @brief Abort with the cuFFT status code if `e` is not `CUFFT_SUCCESS`. */
inline void cufft_ok(cufftResult e, const char *msg) {
  if (e != CUFFT_SUCCESS) {
    std::fprintf(stderr, "cuFFT error: code %d (%s)\n", static_cast<int>(e),
                 msg);
    std::abort();
  }
}

#endif

#endif // DIFFUSION3D_CUDA_CHECK_H
