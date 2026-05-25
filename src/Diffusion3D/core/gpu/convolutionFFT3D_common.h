#ifndef CONVOLUTIONFFT3D_COMMON_H
#define CONVOLUTIONFFT3D_COMMON_H

/**
 * @file convolutionFFT3D_common.h
 * @brief Shared API for 3D FFT-based convolution helpers (GPU).
 *
 * These helpers are based off the 2D convolution helpers in the NVIDIA CUDALibrarySamples.
 *
 * References:
 * - https://github.com/NVIDIA/CUDALibrarySamples
 */

typedef unsigned int uint;

////////////////////////////////////////////////////////////////////////////////
// Helper functions (match CUDA sample style)
////////////////////////////////////////////////////////////////////////////////
// Round a / b to nearest higher integer value
inline int iDivUp(int a, int b) { return (a % b != 0) ? (a / b + 1) : (a / b); }

// Align a to nearest higher multiple of b
inline int iAlignUp(int a, int b) { return (a % b != 0) ? (a - a % b + b) : a; }

#ifdef __CUDACC__
typedef float2 fComplex;
#else
typedef struct
{
    float x;
    float y;
} fComplex;
#endif

extern "C"
{
    /**
     * @brief Pad a small 3D kernel into a larger FFT buffer with ifftshift-on-write.
     *
     * Kernel voxel (z,y,x) is written to:
     *   ((z - cZ + fftZ) % fftZ, (y - cY + fftY) % fftY, (x - cX + fftX) % fftX)
     *
     * @param d_PaddedKernel dst, size fftZ*fftY*fftX (real)
     * @param d_Kernel       src, size kZ*kY*kX (real)
     * @param fftZ,fftY,fftX padded FFT dimensions
     * @param kZ,kY,kX       kernel dimensions
     * @param cZ,cY,cX       kernel center indices (typically k* / 2)
     */
    void padKernel3D(float *d_PaddedKernel,
                     const float *d_Kernel,
                     int fftZ, int fftY, int fftX,
                     int kZ, int kY, int kX,
                     int cZ, int cY, int cX);

    /**
     * @brief Pad 3D data into an FFT buffer using clamp-to-border semantics.
     *
     * @param d_PaddedData dst, size fftZ*fftY*fftX (real)
     * @param d_Data       src, size dZ*dY*dX (real)
     * @param fftZ,fftY,fftX padded FFT dimensions
     * @param dZ,dY,dX       data dimensions
     * @param kZ,kY,kX       kernel dimensions (used to define the border region)
     * @param cZ,cY,cX       kernel center indices
     */
    void padDataClampToBorder3D(float *d_PaddedData,
                                const float *d_Data,
                                int fftZ, int fftY, int fftX,
                                int dZ, int dY, int dX,
                                int kZ, int kY, int kX,
                                int cZ, int cY, int cX);

    /**
     * @brief Extract the physical domain from the padded inverse-transform result.
     *
     * @param d_Dst dst, size dZ*dY*dX (real)
     * @param d_Src src, size fftZ*fftY*fftX (real)
     * @param dZ,dY,dX data dimensions
     * @param fftZ,fftY,fftX padded FFT dimensions
     */
    void unpadResult3D(float *d_Dst,
                       const float *d_Src,
                       int dZ, int dY, int dX,
                       int fftZ, int fftY, int fftX);

    /**
     * @brief Pointwise multiply data spectrum by kernel spectrum and normalize.
     *
     * Mirrors the CUDA sample `modulateAndNormalize` signature:
     * - The packed R2C spectrum has width `(fftX/2 + padding)` along X.
     * - For standard cuFFT R2C output, use `padding = 1` so width is `(fftX/2 + 1)`.
     *
     * Internally, this computes:
     *   elements = fftZ * fftY * (fftX/2 + padding)
     * and applies normalization by `1/(fftZ*fftY*fftX)`.
     */
    void modulateAndNormalize3D(fComplex *d_Dst,
                                fComplex *d_Src,
                                int fftZ,
                                int fftY,
                                int fftX,
                                int padding);

    /**
     * @brief Pointwise integer power of each R2C spectrum bin: dst[i] = src[i]^power (complex).
     *
     * Packed layout matches `modulateAndNormalize3D` with the same padding flag (typically 1).
     * For power 1, copies src to dst. Used to form (K_hat)^N.
     */
    void spectrumRaisePowInt3D(fComplex *d_Dst,
                               const fComplex *d_Src,
                               int fftZ,
                               int fftY,
                               int fftX,
                               int padding,
                               int power);
}

#endif // CONVOLUTIONFFT3D_COMMON_H

