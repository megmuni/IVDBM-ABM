/**
 * @file convolutionFFT3D.cu
 * @brief CUDA-side building blocks for 3D FFT-based convolution helpers.
 *
 * Based on the 2D convolution helpers in the NVIDIA CUDALibrarySamples.
 *
 * References:
 * - https://github.com/NVIDIA/CUDALibrarySamples
 */

#include <assert.h>
#include <stdio.h>

#include "convolutionFFT3D_common.h"

#ifdef __CUDACC__

#include <cuda_runtime.h>

#include "diffusion3d_cuda_check.h"

static inline void cudaCheckLastError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        printf("CUDA error: %s: %s\n", msg, cudaGetErrorString(err));
        assert(false);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Spectrum modulation + normalization.
//
// Thread mapping: 1 thread -> 1 packed complex element in the R2C spectrum.
// Semantics: d_DataSpectrum[i] = d_DataSpectrum[i] * d_KernelSpectrum[i] * scale.
////////////////////////////////////////////////////////////////////////////////
inline __device__ void mulAndScale(fComplex &a, const fComplex &b, const float &c)
{
    fComplex t = {c * (a.x * b.x - a.y * b.y), c * (a.y * b.x + a.x * b.y)};
    a = t;
}

__global__ void modulateAndNormalize3D_kernel(fComplex *d_Dst,
                                              const fComplex *d_Src,
                                              int elements,
                                              float scale)
{
    const int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= elements) return;

    fComplex a = d_Src[i];
    fComplex b = d_Dst[i];
    mulAndScale(a, b, scale);
    d_Dst[i] = a;
}

////////////////////////////////////////////////////////////////////////////////
// Unpad result: extract the physical domain from a padded volume.
//
// Thread mapping: 1 thread -> 1 output voxel in the physical domain.
// This matches the CUDA sample convention used for CPU reference comparisons:
// take the top-left-front [0..dZ)×[0..dY)×[0..dX) window.
////////////////////////////////////////////////////////////////////////////////
__global__ void unpadResult3D_kernel(float *d_Dst,
                                     const float *d_Src,
                                     int dZ,
                                     int dY,
                                     int dX,
                                     int fftZ,
                                     int fftY,
                                     int fftX)
{
    const int x = blockDim.x * blockIdx.x + threadIdx.x;
    const int y = blockDim.y * blockIdx.y + threadIdx.y;
    const int z = blockDim.z * blockIdx.z + threadIdx.z;

    if (x >= dX || y >= dY || z >= dZ) return;

    const int src_idx = (z * fftY + y) * fftX + x;
    const int dst_idx = (z * dY + y) * dX + x;
    d_Dst[dst_idx] = d_Src[src_idx];
}

extern "C" void unpadResult3D(float *d_Dst,
                              const float *d_Src,
                              int dZ, int dY, int dX,
                              int fftZ, int fftY, int fftX)
{
    assert(dZ >= 0 && dY >= 0 && dX >= 0);
    assert(fftZ >= dZ && fftY >= dY && fftX >= dX);
    assert(d_Dst != nullptr);
    assert(d_Src != nullptr);

    if (dZ == 0 || dY == 0 || dX == 0) return;

    dim3 threads(8, 8, 8);
    dim3 grid(iDivUp(dX, threads.x), iDivUp(dY, threads.y), iDivUp(dZ, threads.z));
    unpadResult3D_kernel<<<grid, threads>>>(d_Dst, d_Src, dZ, dY, dX, fftZ, fftY, fftX);
    cudaCheckLastError("unpadResult3D_kernel<<<>>>");
}

////////////////////////////////////////////////////////////////////////////////
// Pad kernel: ifftshift-on-write (kernel center -> padded origin).
//
// Thread mapping: 1 thread -> 1 source kernel voxel.
// Destination index uses modular wrap so that the kernel center (cZ,cY,cX)
// becomes (0,0,0) in the padded volume; this avoids an extra spatial shift after
// the inverse FFT.
////////////////////////////////////////////////////////////////////////////////
__global__ void padKernel3D_kernel(float *d_Dst,
                                   const float *d_Src,
                                   int fftZ,
                                   int fftY,
                                   int fftX,
                                   int kZ,
                                   int kY,
                                   int kX,
                                   int cZ,
                                   int cY,
                                   int cX)
{
    const int x = blockDim.x * blockIdx.x + threadIdx.x;
    const int y = blockDim.y * blockIdx.y + threadIdx.y;
    const int z = blockDim.z * blockIdx.z + threadIdx.z;

    if (x >= kX || y >= kY || z >= kZ) return;

    int pz = z - cZ;
    if (pz < 0) pz += fftZ;
    int py = y - cY;
    if (py < 0) py += fftY;
    int px = x - cX;
    if (px < 0) px += fftX;

    const int src_idx = (z * kY + y) * kX + x;
    const int dst_idx = (pz * fftY + py) * fftX + px;
    d_Dst[dst_idx] = d_Src[src_idx];
}

extern "C" void padKernel3D(float *d_PaddedKernel,
                            const float *d_Kernel,
                            int fftZ, int fftY, int fftX,
                            int kZ, int kY, int kX,
                            int cZ, int cY, int cX)
{
    assert(d_PaddedKernel != nullptr);
    assert(d_Kernel != nullptr);
    assert(fftZ >= kZ && fftY >= kY && fftX >= kX);
    assert(cZ >= 0 && cZ < kZ);
    assert(cY >= 0 && cY < kY);
    assert(cX >= 0 && cX < kX);

    if (kZ == 0 || kY == 0 || kX == 0) return;

    dim3 threads(8, 8, 8);
    dim3 grid(iDivUp(kX, threads.x), iDivUp(kY, threads.y), iDivUp(kZ, threads.z));
    padKernel3D_kernel<<<grid, threads>>>(d_PaddedKernel,
                                          d_Kernel,
                                          fftZ, fftY, fftX,
                                          kZ, kY, kX,
                                          cZ, cY, cX);
    cudaCheckLastError("padKernel3D_kernel<<<>>>");
}

////////////////////////////////////////////////////////////////////////////////
// Pad data with clamp-to-border semantics (CUDA sample behavior).
//
// Thread mapping: 1 thread -> 1 padded voxel.
//
// For each axis independently (shown for z):
//   borderZ = dZ + cZ
//   if z < dZ              -> dz = z
//   else if z < borderZ    -> dz = dZ - 1
//   else                   -> dz = 0
//
// This is not a pure clamp; it clamps only inside the "border region" and then
// wraps to 0 beyond it, matching the 2D CUDA sample's padDataClampToBorder.
////////////////////////////////////////////////////////////////////////////////
__global__ void padDataClampToBorder3D_kernel(float *d_Dst,
                                              const float *d_Src,
                                              int fftZ,
                                              int fftY,
                                              int fftX,
                                              int dZ,
                                              int dY,
                                              int dX,
                                              int cZ,
                                              int cY,
                                              int cX)
{
    const int x = blockDim.x * blockIdx.x + threadIdx.x;
    const int y = blockDim.y * blockIdx.y + threadIdx.y;
    const int z = blockDim.z * blockIdx.z + threadIdx.z;

    if (x >= fftX || y >= fftY || z >= fftZ) return;

    const int borderZ = dZ + cZ;
    const int borderY = dY + cY;
    const int borderX = dX + cX;

    int dz;
    if (z < dZ) dz = z;
    else if (z < borderZ) dz = dZ - 1;
    else dz = 0;

    int dy;
    if (y < dY) dy = y;
    else if (y < borderY) dy = dY - 1;
    else dy = 0;

    int dx;
    if (x < dX) dx = x;
    else if (x < borderX) dx = dX - 1;
    else dx = 0;

    const int src_idx = (dz * dY + dy) * dX + dx;
    const int dst_idx = (z * fftY + y) * fftX + x;
    d_Dst[dst_idx] = d_Src[src_idx];
}

extern "C" void padDataClampToBorder3D(float *d_PaddedData,
                                       const float *d_Data,
                                       int fftZ, int fftY, int fftX,
                                       int dZ, int dY, int dX,
                                       int kZ, int kY, int kX,
                                       int cZ, int cY, int cX)
{
    assert(d_PaddedData != nullptr);
    assert(d_Data != nullptr);
    assert(fftZ >= 0 && fftY >= 0 && fftX >= 0);
    assert(dZ >= 0 && dY >= 0 && dX >= 0);
    assert(fftZ >= dZ && fftY >= dY && fftX >= dX);
    assert(kZ >= 0 && kY >= 0 && kX >= 0);
    assert(cZ >= 0 && cZ < kZ);
    assert(cY >= 0 && cY < kY);
    assert(cX >= 0 && cX < kX);
    // Sample-style border region: requires enough padding to reach (d* + c*).
    assert(fftZ >= dZ + cZ);
    assert(fftY >= dY + cY);
    assert(fftX >= dX + cX);

    if (fftZ == 0 || fftY == 0 || fftX == 0) return;
    if (dZ == 0 || dY == 0 || dX == 0) return;

    dim3 threads(8, 8, 8);
    dim3 grid(iDivUp(fftX, threads.x), iDivUp(fftY, threads.y), iDivUp(fftZ, threads.z));
    padDataClampToBorder3D_kernel<<<grid, threads>>>(d_PaddedData,
                                                     d_Data,
                                                     fftZ, fftY, fftX,
                                                     dZ, dY, dX,
                                                     cZ, cY, cX);
    cudaCheckLastError("padDataClampToBorder3D_kernel<<<>>>");
}

////////////////////////////////////////////////////////////////////////////////
// Pad data by mirroring, for reflecting (zero-flux Neumann) boundaries.
//
// Thread mapping: 1 thread -> 1 padded voxel.
//
// The padded box is exactly 2*d per axis and holds the even reflection of the
// domain about its outer faces (half-sample symmetric):
//
//   dz = z            for z <  dZ
//   dz = 2*dZ - 1 - z for z >= dZ
//
// Under the FFT's periodic wrap this reproduces the infinite train of mirror
// images that the method of images places behind a reflecting wall, so the
// circular convolution equals the reflecting-BC result on the first d cells:
//
//   v[-1]  == v[2*dZ - 1] == u[0]      (ghost below the low face)
//   v[dZ]              == u[dZ - 1]    (ghost above the high face)
//
// which is exactly the ghost value a zero-flux boundary wants. Because the
// extension is the true image train and not a finite skirt, this holds at ANY
// kernel radius -- unlike padDataClampToBorder3D, whose border region is only
// c* voxels wide and is therefore correct only for a radius-1 kernel.
////////////////////////////////////////////////////////////////////////////////
__global__ void padDataMirror3D_kernel(float *d_Dst,
                                       const float *d_Src,
                                       int fftZ,
                                       int fftY,
                                       int fftX,
                                       int dZ,
                                       int dY,
                                       int dX)
{
    const int x = blockDim.x * blockIdx.x + threadIdx.x;
    const int y = blockDim.y * blockIdx.y + threadIdx.y;
    const int z = blockDim.z * blockIdx.z + threadIdx.z;

    if (x >= fftX || y >= fftY || z >= fftZ) return;

    const int dz = (z < dZ) ? z : (2 * dZ - 1 - z);
    const int dy = (y < dY) ? y : (2 * dY - 1 - y);
    const int dx = (x < dX) ? x : (2 * dX - 1 - x);

    const int src_idx = (dz * dY + dy) * dX + dx;
    const int dst_idx = (z * fftY + y) * fftX + x;
    d_Dst[dst_idx] = d_Src[src_idx];
}

extern "C" void padDataMirror3D(float *d_PaddedData,
                                const float *d_Data,
                                int fftZ, int fftY, int fftX,
                                int dZ, int dY, int dX)
{
    assert(d_PaddedData != nullptr);
    assert(d_Data != nullptr);
    assert(dZ >= 0 && dY >= 0 && dX >= 0);
    // The mirror is only a mirror when the box is exactly twice the domain; any
    // other extent would leave a gap that reintroduces a wrong boundary.
    assert(fftZ == 2 * dZ);
    assert(fftY == 2 * dY);
    assert(fftX == 2 * dX);

    if (dZ == 0 || dY == 0 || dX == 0) return;

    dim3 threads(8, 8, 8);
    dim3 grid(iDivUp(fftX, threads.x), iDivUp(fftY, threads.y), iDivUp(fftZ, threads.z));
    padDataMirror3D_kernel<<<grid, threads>>>(d_PaddedData,
                                              d_Data,
                                              fftZ, fftY, fftX,
                                              dZ, dY, dX);
    cudaCheckLastError("padDataMirror3D_kernel<<<>>>");
}

extern "C" void modulateAndNormalize3D(fComplex *d_DataSpectrum,
                                       fComplex *d_KernelSpectrum,
                                       int fftZ,
                                       int fftY,
                                       int fftX,
                                       int padding)
{
    // Mirror the 2D helper's shape assumptions in debug builds.
    assert(fftX % 2 == 0);
    assert(fftZ >= 0 && fftY >= 0 && fftX >= 0);
    assert(d_DataSpectrum != nullptr);
    assert(d_KernelSpectrum != nullptr);

    // Packed R2C spectrum width is (fftX/2 + 1). The sample passes padding=1.
    const int dataSize = fftZ * fftY * (fftX / 2 + padding);
    if (dataSize == 0) return;

    // cuFFT is unnormalized; apply exactly one 1/(N) scale.
    const float scale = 1.0f / (float)(fftZ * fftY * fftX);

    modulateAndNormalize3D_kernel<<<iDivUp(dataSize, 256), 256>>>(
        d_DataSpectrum,
        d_KernelSpectrum,
        dataSize,
        scale);
    cudaCheckLastError("modulateAndNormalize3D_kernel<<<>>>");
}

////////////////////////////////////////////////////////////////////////////////
// Spectrum integer power (per-bin exponentiation by squaring).
//
// Thread mapping: 1 thread -> 1 packed complex element.
//
// Why raising the spectrum is valid:
// - Convolution theorem: applying a convolution in real space multiplies spectra in frequency space.
// - Repeating the same convolution N times multiplies by the same spectrum N times, i.e. H(omega)^N.
// References:
// - https://en.wikipedia.org/wiki/Convolution_theorem
//
// Why exponentiation by squaring:
// - We compute z^N in O(log N) complex multiplies per bin instead of O(N).
// Reference for algorithm:
// - https://en.wikipedia.org/wiki/Exponentiation_by_squaring
////////////////////////////////////////////////////////////////////////////////

__device__ fComplex spectrum_cmul(fComplex a, fComplex b)
{
    fComplex t;
    t.x = a.x * b.x - a.y * b.y;
    t.y = a.y * b.x + a.x * b.y;
    return t;
}

// Algorithm: https://en.wikipedia.org/wiki/Exponentiation_by_squaring
__device__ fComplex spectrum_cpow_u32(fComplex z, unsigned int n)
{
    fComplex r = {1.f, 0.f};
    fComplex b = z;
    while (n > 0)
    {
        if (n & 1u)
            r = spectrum_cmul(r, b);
        b = spectrum_cmul(b, b);
        n >>= 1u;
    }
    return r;
}

// The identity element of the complex number multiplication is 1 + 0i.
__global__ void spectrumOnes3D_kernel(fComplex *d_dst, int elements)
{
    const int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= elements)
        return;
    d_dst[i].x = 1.f;
    d_dst[i].y = 0.f;
}

// Perform exponentiation by squaring for each complex number in the spectrum.
__global__ void spectrumRaisePowInt3D_kernel(fComplex *d_dst,
                                            const fComplex *d_src,
                                            int elements,
                                            unsigned int power)
{
    const int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= elements)
        return;
    d_dst[i] = spectrum_cpow_u32(d_src[i], power);
}

// Handles when power is 0 and 1 separately.
// - When power is 0, set all elements to 1 + 0i.
// - When power is 1, copy the source spectrum to the destination spectrum.
// - For other powers, perform exponentiation by squaring for each complex number in the spectrum.
extern "C" void spectrumRaisePowInt3D(fComplex *d_Dst,
                                      const fComplex *d_Src,
                                      int fftZ,
                                      int fftY,
                                      int fftX,
                                      int padding,
                                      int power)
{
    assert(fftX % 2 == 0);
    assert(d_Dst != nullptr);
    assert(d_Src != nullptr);
    assert(power >= 0);

    const int dataSize = fftZ * fftY * (fftX / 2 + padding);
    if (dataSize == 0)
        return;

    if (power == 0)
    {
        spectrumOnes3D_kernel<<<iDivUp(dataSize, 256), 256>>>(d_Dst, dataSize);
        cudaCheckLastError("spectrumOnes3D_kernel<<<>>>");
        return;
    }

    if (power == 1)
    {
        cuda_ok(cudaMemcpy(d_Dst, d_Src, static_cast<size_t>(dataSize) * sizeof(fComplex),
                           cudaMemcpyDeviceToDevice),
                "spectrumRaisePowInt3D D2D copy (power=1)");
        return;
    }

    spectrumRaisePowInt3D_kernel<<<iDivUp(dataSize, 256), 256>>>(
        d_Dst, d_Src, dataSize, static_cast<unsigned int>(power));
    cudaCheckLastError("spectrumRaisePowInt3D_kernel<<<>>>");
}

#endif // __CUDACC__

