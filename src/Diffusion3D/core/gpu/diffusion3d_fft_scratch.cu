/**
 * @file diffusion3d_fft_scratch.cu
 * @brief Device-backed FFT scratch: padded grid convolution with (K_hat)^N combined spectrum per advance.
 */

#include <cassert>
#include <cstddef>
#include <vector>

#include <cuda_runtime.h>
#include <cufft.h>

#include "convolutionFFT3D_common.h"
#include "diffusion3d_cuda_check.h"
#include "diffusion3d_fft_scratch.h"
#include "diffusion3d_timestep.h"

struct DiffusionFftScratch
{
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int fft_x = 0;
    int fft_y = 0;
    int fft_z = 0;

    cufftHandle plan_r2c{};
    cufftHandle plan_c2r{};

    float *d_data = nullptr;
    float *d_data_pad = nullptr;
    float *d_out = nullptr;
    float *d_kernel = nullptr;
    float *d_kernel_pad = nullptr;

    fComplex *d_data_spec = nullptr;
    fComplex *d_kernel_hat_one_step = nullptr; // R2C of padded explicit-Euler mask for current dt_sub
    fComplex *d_kernel_hat_composed = nullptr; // (d_kernel_hat_one_step)^n_sub pointwise
};

DiffusionFftScratch *diffusion3d_fft_scratch_create(int nx, int ny, int nz,
                                                    int fft_x, int fft_y, int fft_z)
{
    assert(nx > 0 && ny > 0 && nz > 0);
    assert((fft_x % 2) == 0);

    const int kX = 3, kY = 3, kZ = 3;
    const int cX = 1, cY = 1, cZ = 1;
    assert(fft_x >= nx + cX && fft_y >= ny + cY && fft_z >= nz + cZ);

    auto *s = new DiffusionFftScratch{};
    s->nx = nx;
    s->ny = ny;
    s->nz = nz;
    s->fft_x = fft_x;
    s->fft_y = fft_y;
    s->fft_z = fft_z;

    const int fft_n = fft_x * fft_y * fft_z;
    const int spec_x = fft_x / 2 + 1;
    const int spec_n = fft_z * fft_y * spec_x;
    const int n = nx * ny * nz;

    cuda_ok(cudaMalloc(reinterpret_cast<void **>(&s->d_data), static_cast<size_t>(n) * sizeof(float)),
            "cudaMalloc d_data");
    cuda_ok(cudaMalloc(reinterpret_cast<void **>(&s->d_out), static_cast<size_t>(n) * sizeof(float)),
            "cudaMalloc d_out");
    cuda_ok(cudaMalloc(reinterpret_cast<void **>(&s->d_kernel),
                       static_cast<size_t>(kX * kY * kZ) * sizeof(float)),
            "cudaMalloc d_kernel");
    cuda_ok(cudaMalloc(reinterpret_cast<void **>(&s->d_data_pad),
                       static_cast<size_t>(fft_n) * sizeof(float)),
            "cudaMalloc d_data_pad");
    cuda_ok(cudaMalloc(reinterpret_cast<void **>(&s->d_kernel_pad),
                       static_cast<size_t>(fft_n) * sizeof(float)),
            "cudaMalloc d_kernel_pad");
    cuda_ok(cudaMalloc(reinterpret_cast<void **>(&s->d_data_spec),
                       static_cast<size_t>(spec_n) * sizeof(fComplex)),
            "cudaMalloc d_data_spec");
    cuda_ok(cudaMalloc(reinterpret_cast<void **>(&s->d_kernel_hat_one_step),
                       static_cast<size_t>(spec_n) * sizeof(fComplex)),
            "cudaMalloc d_kernel_hat_one_step");
    cuda_ok(cudaMalloc(reinterpret_cast<void **>(&s->d_kernel_hat_composed),
                       static_cast<size_t>(spec_n) * sizeof(fComplex)),
            "cudaMalloc d_kernel_hat_composed");

    cufft_ok(cufftPlan3d(&s->plan_r2c, fft_z, fft_y, fft_x, CUFFT_R2C), "cufftPlan3d R2C");
    cufft_ok(cufftPlan3d(&s->plan_c2r, fft_z, fft_y, fft_x, CUFFT_C2R), "cufftPlan3d C2R");

    return s;
}

void diffusion3d_fft_scratch_destroy(DiffusionFftScratch *s)
{
    if (s == nullptr)
        return;

    cufftDestroy(s->plan_r2c);
    cufftDestroy(s->plan_c2r);

    cudaFree(s->d_data);
    cudaFree(s->d_out);
    cudaFree(s->d_kernel);
    cudaFree(s->d_data_pad);
    cudaFree(s->d_kernel_pad);
    cudaFree(s->d_data_spec);
    cudaFree(s->d_kernel_hat_one_step);
    cudaFree(s->d_kernel_hat_composed);

    delete s;
}

SubstepPlan diffusion3d_fft_scratch_rebuild_operator(DiffusionFftScratch *s,
                                                     double h,
                                                     double D,
                                                     double tick_dt,
                                                     double safety)
{
    assert(s != nullptr);
    assert(tick_dt >= 0.0);
    assert(safety > 0.0 && safety <= 1.0);

    const int fft_x = s->fft_x;
    const int fft_y = s->fft_y;
    const int fft_z = s->fft_z;

    const double dt_max = compute_stability_constraint(h, D, safety);
    const SubstepPlan plan = plan_substeps(tick_dt, dt_max);

    const int kX = 3, kY = 3, kZ = 3;
    const int cX = 1, cY = 1, cZ = 1;

    const double dt_sub = plan.dt_sub;
    const double s_fac = dt_sub * D / (h * h);

    std::vector<float> h_kernel(static_cast<size_t>(kX * kY * kZ), 0.0f);
    auto kidx = [&](int z, int y, int x) { return (z * kY + y) * kX + x; };
    h_kernel[static_cast<size_t>(kidx(cZ, cY, cX))] = static_cast<float>(1.0 - 6.0 * s_fac);
    h_kernel[static_cast<size_t>(kidx(cZ, cY, cX + 1))] = static_cast<float>(s_fac);
    h_kernel[static_cast<size_t>(kidx(cZ, cY, cX - 1))] = static_cast<float>(s_fac);
    h_kernel[static_cast<size_t>(kidx(cZ, cY + 1, cX))] = static_cast<float>(s_fac);
    h_kernel[static_cast<size_t>(kidx(cZ, cY - 1, cX))] = static_cast<float>(s_fac);
    h_kernel[static_cast<size_t>(kidx(cZ + 1, cY, cX))] = static_cast<float>(s_fac);
    h_kernel[static_cast<size_t>(kidx(cZ - 1, cY, cX))] = static_cast<float>(s_fac);

    cuda_ok(cudaMemcpy(s->d_kernel, h_kernel.data(), h_kernel.size() * sizeof(float),
                       cudaMemcpyHostToDevice),
            "H2D kernel");

    const int fft_n = fft_x * fft_y * fft_z;
    cuda_ok(cudaMemset(s->d_kernel_pad, 0, static_cast<size_t>(fft_n) * sizeof(float)), "memset kernel_pad");
    padKernel3D(s->d_kernel_pad, s->d_kernel,
                fft_z, fft_y, fft_x,
                kZ, kY, kX,
                cZ, cY, cX);
    cudaDeviceSynchronize();

    cufft_ok(cufftExecR2C(s->plan_r2c, reinterpret_cast<cufftReal *>(s->d_kernel_pad),
                          reinterpret_cast<cufftComplex *>(s->d_kernel_hat_one_step)),
             "cufftExecR2C kernel_pad");
    cudaDeviceSynchronize();

    spectrumRaisePowInt3D(s->d_kernel_hat_composed, s->d_kernel_hat_one_step,
                          fft_z, fft_y, fft_x,
                          /*padding=*/1,
                          plan.n_sub);
    cudaDeviceSynchronize();

    return plan;
}

void diffusion3d_fft_scratch_apply(DiffusionFftScratch *s, std::vector<double> &field)
{
    assert(s != nullptr);
    assert(static_cast<int>(field.size()) == s->nx * s->ny * s->nz);

    const int fft_x = s->fft_x;
    const int fft_y = s->fft_y;
    const int fft_z = s->fft_z;
    const int fft_n = fft_x * fft_y * fft_z;
    const int n = s->nx * s->ny * s->nz;

    std::vector<float> h_tmp(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        h_tmp[static_cast<size_t>(i)] = static_cast<float>(field[static_cast<size_t>(i)]);

    cuda_ok(cudaMemcpy(s->d_data, h_tmp.data(), static_cast<size_t>(n) * sizeof(float),
                       cudaMemcpyHostToDevice),
            "H2D field");

    cuda_ok(cudaMemset(s->d_data_pad, 0, static_cast<size_t>(fft_n) * sizeof(float)), "memset data_pad");
    // Mirror, not clamp-to-border: the composed operator has radius n_sub, but the
    // clamp's border region is only one voxel wide, so every ghost layer past the
    // first read slice 0 and corrupted the boundary. The mirror is the method-of-
    // images extension and is exact at any radius.
    padDataMirror3D(s->d_data_pad, s->d_data,
                    fft_z, fft_y, fft_x,
                    s->nz, s->ny, s->nx);
    cudaDeviceSynchronize();

    cufft_ok(cufftExecR2C(s->plan_r2c, reinterpret_cast<cufftReal *>(s->d_data_pad),
                          reinterpret_cast<cufftComplex *>(s->d_data_spec)),
             "cufftExecR2C data_pad");
    cudaDeviceSynchronize();

    modulateAndNormalize3D(s->d_data_spec, s->d_kernel_hat_composed,
                           fft_z, fft_y, fft_x,
                           /*padding=*/1);
    cudaDeviceSynchronize();

    cufft_ok(cufftExecC2R(s->plan_c2r, reinterpret_cast<cufftComplex *>(s->d_data_spec),
                          reinterpret_cast<cufftReal *>(s->d_data_pad)),
             "cufftExecC2R data_spec");
    cudaDeviceSynchronize();

    unpadResult3D(s->d_out, s->d_data_pad, s->nz, s->ny, s->nx, fft_z, fft_y, fft_x);
    cudaDeviceSynchronize();

    cuda_ok(cudaMemcpy(h_tmp.data(), s->d_out, static_cast<size_t>(n) * sizeof(float),
                       cudaMemcpyDeviceToHost),
            "D2H field");

    for (int i = 0; i < n; ++i)
    {
        const double v = static_cast<double>(h_tmp[static_cast<size_t>(i)]);
        field[static_cast<size_t>(i)] = (v > 0.0) ? v : 0.0;
    }
}

