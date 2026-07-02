/**
 * @file convolutionFFT3D_test.cpp
 * @brief Catch2 tests for `modulateAndNormalize3D`.
 */

#include <vector>
#include <cmath>

#include <catch2/catch.hpp>

#ifdef DIFFUSION3D_CUDA
#include <cuda_runtime.h>
#include <cufft.h>
#endif

#include "convolutionFFT3D_common.h"
#include "analytic_diffusion3d.h"
#include "diffusion3d_step_euler_cpu.h"
#include "scalar_field_grid.h"

static inline bool nearly_equal(float a, float b, float eps = 1e-6f)
{
    return std::fabs(a - b) <= eps;
}

TEST_CASE("modulateAndNormalize3D matches expected output", "[fft][helpers][modulate3d]")
{
#ifndef DIFFUSION3D_CUDA
    SUCCEED("DIFFUSION3D_CUDA is OFF; skipping CUDA helper test.");
#else
    // Mirror the 2D sample interface: elements = fftZ*fftY*(fftX/2 + padding)
    const int fftZ = 1;
    const int fftY = 1;
    const int fftX = 32;
    const int padding = 1;
    const int elements = fftZ * fftY * (fftX / 2 + padding);
    const float scale = 1.0f / float(fftZ * fftY * fftX);

    std::vector<fComplex> h_dst(elements);
    std::vector<fComplex> h_src(elements);
    std::vector<fComplex> h_expected(elements);

    for (int i = 0; i < elements; i++)
    {
        h_dst[i].x = 0.1f * (i + 1);
        h_dst[i].y = -0.05f * (i + 2);
        h_src[i].x = -0.2f * (i + 3);
        h_src[i].y = 0.07f * (i + 4);

        // Expected: (src * dst) * scale
        const float ax = h_src[i].x, ay = h_src[i].y;
        const float bx = h_dst[i].x, by = h_dst[i].y;
        h_expected[i].x = scale * (ax * bx - ay * by);
        h_expected[i].y = scale * (ay * bx + ax * by);
    }

    fComplex *d_dst = nullptr;
    fComplex *d_src = nullptr;
    cudaMalloc(&d_dst, sizeof(fComplex) * elements);
    cudaMalloc(&d_src, sizeof(fComplex) * elements);

    cudaMemcpy(d_dst, h_dst.data(), sizeof(fComplex) * elements, cudaMemcpyHostToDevice);
    cudaMemcpy(d_src, h_src.data(), sizeof(fComplex) * elements, cudaMemcpyHostToDevice);

    modulateAndNormalize3D(d_dst, d_src, fftZ, fftY, fftX, padding);
    cudaDeviceSynchronize();

    std::vector<fComplex> h_out(elements);
    cudaMemcpy(h_out.data(), d_dst, sizeof(fComplex) * elements, cudaMemcpyDeviceToHost);

    cudaFree(d_dst);
    cudaFree(d_src);

    for (int i = 0; i < elements; i++)
    {
        INFO("i=" << i
                 << " out=(" << h_out[i].x << "," << h_out[i].y << ")"
                 << " expected=(" << h_expected[i].x << "," << h_expected[i].y << ")");
        CHECK(nearly_equal(h_out[i].x, h_expected[i].x, 1e-6f));
        CHECK(nearly_equal(h_out[i].y, h_expected[i].y, 1e-6f));
    }
#endif
}

namespace
{

fComplex spectrum_cmul_host(fComplex a, fComplex b)
{
    fComplex t;
    t.x = a.x * b.x - a.y * b.y;
    t.y = a.y * b.x + a.x * b.y;
    return t;
}

fComplex spectrum_cpow_u32_host(fComplex z, unsigned int n)
{
    fComplex r{1.f, 0.f};
    fComplex b = z;
    while (n > 0)
    {
        if (n & 1u)
            r = spectrum_cmul_host(r, b);
        b = spectrum_cmul_host(b, b);
        n >>= 1u;
    }
    return r;
}

} // namespace

TEST_CASE("spectrumRaisePowInt3D matches host complex powers", "[fft][helpers][spectrum_pow]")
{
#ifndef DIFFUSION3D_CUDA
    SUCCEED("DIFFUSION3D_CUDA is OFF; skipping spectrumRaisePowInt3D test.");
#else
    const int fftZ = 2;
    const int fftY = 2;
    const int fftX = 8;
    const int padding = 1;
    const int elements = fftZ * fftY * (fftX / 2 + padding);
    REQUIRE((fftX % 2) == 0);

    std::vector<fComplex> h_src(elements);
    for (int i = 0; i < elements; ++i)
    {
        h_src[static_cast<size_t>(i)].x = 0.02f * static_cast<float>(i + 1) - 0.3f;
        h_src[static_cast<size_t>(i)].y = 0.015f * static_cast<float>(i + 2) + 0.1f;
    }

    fComplex *d_src = nullptr;
    fComplex *d_dst = nullptr;
    REQUIRE(cudaMalloc(reinterpret_cast<void **>(&d_src), sizeof(fComplex) * static_cast<size_t>(elements)) ==
            cudaSuccess);
    REQUIRE(cudaMalloc(reinterpret_cast<void **>(&d_dst), sizeof(fComplex) * static_cast<size_t>(elements)) ==
            cudaSuccess);
    REQUIRE(cudaMemcpy(d_src, h_src.data(), sizeof(fComplex) * static_cast<size_t>(elements),
                       cudaMemcpyHostToDevice) == cudaSuccess);

    for (int power : {0, 1, 2, 3, 7})
    {
        spectrumRaisePowInt3D(d_dst, d_src, fftZ, fftY, fftX, padding, power);
        cudaDeviceSynchronize();

        std::vector<fComplex> h_out(static_cast<size_t>(elements));
        REQUIRE(cudaMemcpy(h_out.data(), d_dst, sizeof(fComplex) * h_out.size(), cudaMemcpyDeviceToHost) ==
                cudaSuccess);

        const float eps = 2e-5f;
        for (int i = 0; i < elements; ++i)
        {
            const fComplex z = h_src[static_cast<size_t>(i)];
            const fComplex e = spectrum_cpow_u32_host(z, static_cast<unsigned int>(power));
            INFO("power=" << power << " i=" << i);
            CHECK(std::fabs(h_out[static_cast<size_t>(i)].x - e.x) < eps);
            CHECK(std::fabs(h_out[static_cast<size_t>(i)].y - e.y) < eps);
        }
    }

    cudaFree(d_dst);
    cudaFree(d_src);
#endif
}

TEST_CASE("Self-convolution: repeated one-step FFT matches spectrum power", "[fft][e2e][cufft][self_convolution]")
{
#ifndef DIFFUSION3D_CUDA
    SUCCEED("DIFFUSION3D_CUDA is OFF; skipping self-convolution test.");
#else
    // This test validates the algebra behind the FFT backend:
    // - Repeating the same one-step convolution N times corresponds to raising the kernel spectrum to N.
    // - We compare the full end-to-end pipeline (pad -> R2C -> modulate -> C2R -> unpad) in both forms.

    const int dX = 6, dY = 6, dZ = 6;
    const int fftX = 16, fftY = 16, fftZ = 16;
    const int kX = 3, kY = 3, kZ = 3;
    const int cX = 1, cY = 1, cZ = 1;
    const int N = 7; // Number of one-step convolutions.

    REQUIRE((fftX % 2) == 0);
    REQUIRE(fftX >= dX + cX);
    REQUIRE(fftY >= dY + cY);
    REQUIRE(fftZ >= dZ + cZ);

    const int dN = dX * dY * dZ;
    const int fftN = fftX * fftY * fftZ;
    const int specX = fftX / 2 + 1;
    const int specN = fftZ * fftY * specX;

    // Use a small explicit-Euler step so kernel stays well-conditioned.
    const double dt_sub = 1e-3;
    const double h = 1.0;
    const double D = 0.1;
    const double s = dt_sub * D / (h * h);

    std::vector<float> h_kernel(static_cast<size_t>(kX * kY * kZ), 0.0f);
    auto kidx = [&](int z, int y, int x) { return (z * kY + y) * kX + x; };
    h_kernel[static_cast<size_t>(kidx(cZ, cY, cX))] = static_cast<float>(1.0 - 6.0 * s);
    h_kernel[static_cast<size_t>(kidx(cZ, cY, cX + 1))] = static_cast<float>(s);
    h_kernel[static_cast<size_t>(kidx(cZ, cY, cX - 1))] = static_cast<float>(s);
    h_kernel[static_cast<size_t>(kidx(cZ, cY + 1, cX))] = static_cast<float>(s);
    h_kernel[static_cast<size_t>(kidx(cZ, cY - 1, cX))] = static_cast<float>(s);
    h_kernel[static_cast<size_t>(kidx(cZ + 1, cY, cX))] = static_cast<float>(s);
    h_kernel[static_cast<size_t>(kidx(cZ - 1, cY, cX))] = static_cast<float>(s);

    // Deterministic input on the physical domain.
    std::vector<float> h_u0(static_cast<size_t>(dN), 0.0f);
    for (int z = 0; z < dZ; ++z)
    {
        for (int y = 0; y < dY; ++y)
        {
            for (int x = 0; x < dX; ++x)
            {
                const int i = (z * dY + y) * dX + x;
                h_u0[static_cast<size_t>(i)] =
                    0.02f * static_cast<float>(i + 1) +
                    0.1f * std::sin(0.3f * static_cast<float>(x)) +
                    0.05f * std::cos(0.2f * static_cast<float>(y + z));
            }
        }
    }

    float *d_data = nullptr;
    float *d_kernel = nullptr;
    float *d_data_pad = nullptr;
    float *d_kernel_pad = nullptr;
    float *d_out = nullptr;
    fComplex *d_data_spec = nullptr;
    fComplex *d_kernel_hat = nullptr;
    fComplex *d_kernel_hat_pow = nullptr;

    REQUIRE(cudaMalloc(&d_data, sizeof(float) * static_cast<size_t>(dN)) == cudaSuccess);
    REQUIRE(cudaMalloc(&d_kernel, sizeof(float) * h_kernel.size()) == cudaSuccess);
    REQUIRE(cudaMalloc(&d_data_pad, sizeof(float) * static_cast<size_t>(fftN)) == cudaSuccess);
    REQUIRE(cudaMalloc(&d_kernel_pad, sizeof(float) * static_cast<size_t>(fftN)) == cudaSuccess);
    REQUIRE(cudaMalloc(&d_out, sizeof(float) * static_cast<size_t>(dN)) == cudaSuccess);
    REQUIRE(cudaMalloc(&d_data_spec, sizeof(fComplex) * static_cast<size_t>(specN)) == cudaSuccess);
    REQUIRE(cudaMalloc(&d_kernel_hat, sizeof(fComplex) * static_cast<size_t>(specN)) == cudaSuccess);
    REQUIRE(cudaMalloc(&d_kernel_hat_pow, sizeof(fComplex) * static_cast<size_t>(specN)) == cudaSuccess);

    REQUIRE(cudaMemcpy(d_kernel, h_kernel.data(), sizeof(float) * h_kernel.size(), cudaMemcpyHostToDevice) ==
            cudaSuccess);
    REQUIRE(cudaMemset(d_kernel_pad, 0, sizeof(float) * static_cast<size_t>(fftN)) == cudaSuccess);
    padKernel3D(d_kernel_pad, d_kernel, fftZ, fftY, fftX, kZ, kY, kX, cZ, cY, cX);
    cudaDeviceSynchronize();

    cufftHandle planR2C{};
    cufftHandle planC2R{};
    REQUIRE(cufftPlan3d(&planR2C, fftZ, fftY, fftX, CUFFT_R2C) == CUFFT_SUCCESS);
    REQUIRE(cufftPlan3d(&planC2R, fftZ, fftY, fftX, CUFFT_C2R) == CUFFT_SUCCESS);

    REQUIRE(cufftExecR2C(planR2C, reinterpret_cast<cufftReal *>(d_kernel_pad),
                         reinterpret_cast<cufftComplex *>(d_kernel_hat)) == CUFFT_SUCCESS);
    cudaDeviceSynchronize();

    spectrumRaisePowInt3D(d_kernel_hat_pow, d_kernel_hat, fftZ, fftY, fftX, /*padding=*/1, /*power=*/N);
    cudaDeviceSynchronize();

    auto apply_kernel_hat_once = [&](const fComplex *d_hat, const std::vector<float> &h_in, std::vector<float> &h_out)
    {
        REQUIRE(cudaMemcpy(d_data, h_in.data(), sizeof(float) * h_in.size(), cudaMemcpyHostToDevice) == cudaSuccess);
        REQUIRE(cudaMemset(d_data_pad, 0, sizeof(float) * static_cast<size_t>(fftN)) == cudaSuccess);

        padDataClampToBorder3D(d_data_pad, d_data, fftZ, fftY, fftX, dZ, dY, dX, kZ, kY, kX, cZ, cY, cX);
        cudaDeviceSynchronize();

        REQUIRE(cufftExecR2C(planR2C, reinterpret_cast<cufftReal *>(d_data_pad),
                             reinterpret_cast<cufftComplex *>(d_data_spec)) == CUFFT_SUCCESS);
        modulateAndNormalize3D(d_data_spec, const_cast<fComplex *>(d_hat), fftZ, fftY, fftX, /*padding=*/1);
        cudaDeviceSynchronize();

        REQUIRE(cufftExecC2R(planC2R, reinterpret_cast<cufftComplex *>(d_data_spec),
                             reinterpret_cast<cufftReal *>(d_data_pad)) == CUFFT_SUCCESS);
        cudaDeviceSynchronize();

        unpadResult3D(d_out, d_data_pad, dZ, dY, dX, fftZ, fftY, fftX);
        cudaDeviceSynchronize();

        h_out.resize(h_in.size());
        REQUIRE(cudaMemcpy(h_out.data(), d_out, sizeof(float) * h_out.size(), cudaMemcpyDeviceToHost) == cudaSuccess);
    };

    // Reference: chain N one-step convolutions using K_hat.
    std::vector<float> h_chain = h_u0;
    std::vector<float> h_tmp;
    for (int i = 0; i < N; ++i)
    {
        apply_kernel_hat_once(d_kernel_hat, h_chain, h_tmp);
        h_chain.swap(h_tmp);
    }

    // Candidate: single convolution using (K_hat)^N.
    std::vector<float> h_pow;
    apply_kernel_hat_once(d_kernel_hat_pow, h_u0, h_pow);

    // Compare.
    double max_abs = 0.0;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < dN; ++i)
    {
        const double ref = static_cast<double>(h_chain[static_cast<size_t>(i)]);
        const double got = static_cast<double>(h_pow[static_cast<size_t>(i)]);
        const double diff = got - ref;
        max_abs = std::max(max_abs, std::fabs(diff));
        num += diff * diff;
        den += ref * ref;
    }
    const double rel_l2 = (den > 0.0) ? std::sqrt(num / den) : std::sqrt(num);
    INFO("rel_l2=" << rel_l2 << " max_abs=" << max_abs);
    CHECK(rel_l2 < 1e-4);
    CHECK(max_abs < 1e-3);

    cufftDestroy(planR2C);
    cufftDestroy(planC2R);
    cudaFree(d_data);
    cudaFree(d_kernel);
    cudaFree(d_data_pad);
    cudaFree(d_kernel_pad);
    cudaFree(d_out);
    cudaFree(d_data_spec);
    cudaFree(d_kernel_hat);
    cudaFree(d_kernel_hat_pow);
#endif
}

TEST_CASE("CPU diffusion is qualitatively close to analytic Gaussian (early time)", "[diffusion][cpu][analytic]")
{
    const int nx = 33, ny = 33, nz = 33;
    const double dt = 1e-3;
    const double h = 1.0;
    const double D = 0.1;

    ScalarFieldGrid grid(nx, ny, nz);
    const int cx = nx / 2, cy = ny / 2, cz = nz / 2;
    grid.at(cx, cy, cz) = 1.0;

    diffusion3d_step_euler_scalar(grid, h, D, dt);

    const double t = dt;
    std::vector<double> ref = analytic_gaussian_impulse_3d(nx, ny, nz, h, D, t, cx, cy, cz, 1.0);

    const double err = rel_l2_error(grid.data_, ref);
    INFO("rel_l2_error=" << err);
    CHECK(err < 5.0);
}

TEST_CASE("unpadResult3D extracts the top-left-front volume", "[fft][helpers][unpad3d]")
{
#ifndef DIFFUSION3D_CUDA
    SUCCEED("DIFFUSION3D_CUDA is OFF; skipping CUDA helper test.");
#else
    const int fftX = 6, fftY = 5, fftZ = 4;
    const int dX = 3, dY = 2, dZ = 2;

    const int fftN = fftX * fftY * fftZ;
    const int dN = dX * dY * dZ;

    std::vector<float> h_src(fftN, 0.0f);
    for (int z = 0; z < fftZ; z++)
        for (int y = 0; y < fftY; y++)
            for (int x = 0; x < fftX; x++)
                h_src[(z * fftY + y) * fftX + x] = float((z + 1) * 100 + (y + 1) * 10 + (x + 1));

    float *d_src = nullptr;
    float *d_dst = nullptr;
    cudaMalloc(&d_src, sizeof(float) * fftN);
    cudaMalloc(&d_dst, sizeof(float) * dN);
    cudaMemcpy(d_src, h_src.data(), sizeof(float) * fftN, cudaMemcpyHostToDevice);
    cudaMemset(d_dst, 0, sizeof(float) * dN);

    unpadResult3D(d_dst, d_src, dZ, dY, dX, fftZ, fftY, fftX);
    cudaDeviceSynchronize();

    std::vector<float> h_out(dN, 0.0f);
    cudaMemcpy(h_out.data(), d_dst, sizeof(float) * dN, cudaMemcpyDeviceToHost);

    cudaFree(d_src);
    cudaFree(d_dst);

    for (int z = 0; z < dZ; z++)
        for (int y = 0; y < dY; y++)
            for (int x = 0; x < dX; x++)
            {
                const float expected = h_src[(z * fftY + y) * fftX + x];
                const float got = h_out[(z * dY + y) * dX + x];
                INFO("x=" << x << " y=" << y << " z=" << z << " got=" << got << " expected=" << expected);
                CHECK(got == Approx(expected));
            }
#endif
}

TEST_CASE("padKernel3D performs ifftshift-on-write", "[fft][helpers][padKernel3d]")
{
#ifndef DIFFUSION3D_CUDA
    SUCCEED("DIFFUSION3D_CUDA is OFF; skipping CUDA helper test.");
#else
    const int fftX = 5, fftY = 5, fftZ = 5;
    const int kX = 3, kY = 3, kZ = 3;
    const int cX = 1, cY = 1, cZ = 1;

    std::vector<float> h_kernel(kX * kY * kZ, 0.0f);
    // Put a single 1.0 at the kernel center.
    h_kernel[(cZ * kY + cY) * kX + cX] = 1.0f;

    const int fftN = fftX * fftY * fftZ;
    std::vector<float> h_padded(fftN, 0.0f);

    float *d_kernel = nullptr;
    float *d_padded = nullptr;
    cudaMalloc(&d_kernel, sizeof(float) * h_kernel.size());
    cudaMalloc(&d_padded, sizeof(float) * fftN);
    cudaMemcpy(d_kernel, h_kernel.data(), sizeof(float) * h_kernel.size(), cudaMemcpyHostToDevice);
    cudaMemset(d_padded, 0, sizeof(float) * fftN);

    padKernel3D(d_padded, d_kernel, fftZ, fftY, fftX, kZ, kY, kX, cZ, cY, cX);
    cudaDeviceSynchronize();

    cudaMemcpy(h_padded.data(), d_padded, sizeof(float) * fftN, cudaMemcpyDeviceToHost);
    cudaFree(d_kernel);
    cudaFree(d_padded);

    // Center voxel should map to origin (0,0,0).
    CHECK(h_padded[0] == Approx(1.0f));

    // All other values should remain zero.
    int nonzero = 0;
    for (int i = 0; i < fftN; i++)
        if (std::fabs(h_padded[i]) > 1e-7f) nonzero++;
    CHECK(nonzero == 1);
#endif
}

TEST_CASE("padDataClampToBorder3D matches CUDA-sample border semantics", "[fft][helpers][padDataClamp3d]")
{
#ifndef DIFFUSION3D_CUDA
    SUCCEED("DIFFUSION3D_CUDA is OFF; skipping CUDA helper test.");
#else
    const int dX = 2, dY = 2, dZ = 2;
    const int fftX = 4, fftY = 4, fftZ = 4;
    const int kX = 3, kY = 3, kZ = 3;
    const int cX = 1, cY = 1, cZ = 1; // border* = d* + c*

    std::vector<float> h_data(dX * dY * dZ, 0.0f);
    for (int z = 0; z < dZ; z++)
        for (int y = 0; y < dY; y++)
            for (int x = 0; x < dX; x++)
                h_data[(z * dY + y) * dX + x] = float((z + 1) * 100 + (y + 1) * 10 + (x + 1));

    const int fftN = fftX * fftY * fftZ;
    std::vector<float> h_padded(fftN, -1.0f);

    float *d_data = nullptr;
    float *d_padded = nullptr;
    cudaMalloc(&d_data, sizeof(float) * h_data.size());
    cudaMalloc(&d_padded, sizeof(float) * fftN);
    cudaMemcpy(d_data, h_data.data(), sizeof(float) * h_data.size(), cudaMemcpyHostToDevice);
    cudaMemset(d_padded, 0, sizeof(float) * fftN);

    padDataClampToBorder3D(d_padded, d_data,
                           fftZ, fftY, fftX,
                           dZ, dY, dX,
                           kZ, kY, kX,
                           cZ, cY, cX);
    cudaDeviceSynchronize();

    cudaMemcpy(h_padded.data(), d_padded, sizeof(float) * fftN, cudaMemcpyDeviceToHost);
    cudaFree(d_data);
    cudaFree(d_padded);

    auto at = [&](int z, int y, int x) -> float {
        return h_padded[(z * fftY + y) * fftX + x];
    };
    auto src = [&](int z, int y, int x) -> float {
        return h_data[(z * dY + y) * dX + x];
    };

    // Direct region: z,y,x < d*  -> identity mapping.
    CHECK(at(0, 0, 0) == Approx(src(0, 0, 0)));
    CHECK(at(1, 1, 1) == Approx(src(1, 1, 1)));

    // Clamp region: index == d* (but < d*+c*) -> clamp to d*-1.
    CHECK(at(2, 0, 0) == Approx(src(1, 0, 0))); // z clamped
    CHECK(at(0, 2, 1) == Approx(src(0, 1, 1))); // y clamped
    CHECK(at(1, 1, 2) == Approx(src(1, 1, 1))); // x clamped

    // Wrap region: index >= d*+c* -> wraps to 0.
    CHECK(at(3, 0, 0) == Approx(src(0, 0, 0)));
    CHECK(at(0, 3, 0) == Approx(src(0, 0, 0)));
    CHECK(at(0, 0, 3) == Approx(src(0, 0, 0)));
#endif
}

TEST_CASE("3D FFT convolution matches CPU 6-point Laplacian update (interior voxels)", "[fft][e2e][cufft][laplacian]")
{
#ifndef DIFFUSION3D_CUDA
    SUCCEED("DIFFUSION3D_CUDA is OFF; skipping cuFFT end-to-end test.");
#else
    // Compare against one CPU explicit Euler step (`diffusion3d_step_euler_scalar`).
    // We only compare interior voxels so boundary-condition differences do not dominate.
    const int nx = 6, ny = 6, nz = 6;
    const double dt = 1e-3;
    const double h = 1.0;
    const double D = 0.1;
    // Matches explicit Euler: u += dt * (D/h^2) * sum_face (u_nei - u).
    const double s = dt * D / (h * h);

    // Set only one step for diffusion
    ScalarFieldGrid grid(nx, ny, nz);
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
                grid.at(x, y, z) = 0.01 * double((z + 1) * 100 + (y + 1) * 10 + (x + 1));

    const std::vector<double> initial_field = grid.data_;
    diffusion3d_step_euler_scalar(grid, h, D, dt);
    const std::vector<double> cpu_out = grid.data_;

    // Build a 3x3x3 convolution kernel that applies one explicit-Euler update:
    // u_next = u + s * (sum6(nei) - 6*u)
    const int kX = 3, kY = 3, kZ = 3;
    const int cX = 1, cY = 1, cZ = 1;
    std::vector<float> h_kernel(kX * kY * kZ, 0.0f);
    auto kidx = [&](int z, int y, int x) { return (z * kY + y) * kX + x; };
    h_kernel[kidx(cZ, cY, cX)] = float(1.0 - 6.0 * s);
    h_kernel[kidx(cZ, cY, cX + 1)] = float(s);
    h_kernel[kidx(cZ, cY, cX - 1)] = float(s);
    h_kernel[kidx(cZ, cY + 1, cX)] = float(s);
    h_kernel[kidx(cZ, cY - 1, cX)] = float(s);
    h_kernel[kidx(cZ + 1, cY, cX)] = float(s);
    h_kernel[kidx(cZ - 1, cY, cX)] = float(s);

    // FFT sizes: must be >= d* + c* for sample-style padding, and fftX even for R2C packing.
    const int dX = nx, dY = ny, dZ = nz;
    const int fftX = 8, fftY = 8, fftZ = 8;
    REQUIRE((fftX % 2) == 0);
    REQUIRE(fftX >= dX + cX);
    REQUIRE(fftY >= dY + cY);
    REQUIRE(fftZ >= dZ + cZ);

    const int fftN = fftX * fftY * fftZ;
    const int specX = fftX / 2 + 1;
    const int specN = fftZ * fftY * specX;

    std::vector<float> h_data(dX * dY * dZ, 0.0f);
    for (int i = 0; i < (int)h_data.size(); i++)
        h_data[i] = static_cast<float>(initial_field[static_cast<size_t>(i)]);

    float *d_data = nullptr;
    float *d_kernel = nullptr;
    float *d_data_pad = nullptr;
    float *d_kernel_pad = nullptr;
    fComplex *d_data_spec = nullptr;
    fComplex *d_kernel_spec = nullptr;
    float *d_out = nullptr;

    cudaMalloc(&d_data, sizeof(float) * h_data.size());
    cudaMalloc(&d_kernel, sizeof(float) * h_kernel.size());
    cudaMalloc(&d_data_pad, sizeof(float) * fftN);
    cudaMalloc(&d_kernel_pad, sizeof(float) * fftN);
    cudaMalloc(&d_data_spec, sizeof(fComplex) * specN);
    cudaMalloc(&d_kernel_spec, sizeof(fComplex) * specN);
    cudaMalloc(&d_out, sizeof(float) * h_data.size());

    cudaMemcpy(d_data, h_data.data(), sizeof(float) * h_data.size(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_kernel, h_kernel.data(), sizeof(float) * h_kernel.size(), cudaMemcpyHostToDevice);
    cudaMemset(d_data_pad, 0, sizeof(float) * fftN);
    cudaMemset(d_kernel_pad, 0, sizeof(float) * fftN);

    padKernel3D(d_kernel_pad, d_kernel, fftZ, fftY, fftX, kZ, kY, kX, cZ, cY, cX);
    padDataClampToBorder3D(d_data_pad, d_data, fftZ, fftY, fftX, dZ, dY, dX, kZ, kY, kX, cZ, cY, cX);
    cudaDeviceSynchronize();

    cufftHandle planR2C;
    cufftHandle planC2R;
    REQUIRE(cufftPlan3d(&planR2C, fftZ, fftY, fftX, CUFFT_R2C) == CUFFT_SUCCESS);
    REQUIRE(cufftPlan3d(&planC2R, fftZ, fftY, fftX, CUFFT_C2R) == CUFFT_SUCCESS);

    REQUIRE(cufftExecR2C(planR2C, (cufftReal *)d_kernel_pad, (cufftComplex *)d_kernel_spec) == CUFFT_SUCCESS);
    REQUIRE(cufftExecR2C(planR2C, (cufftReal *)d_data_pad, (cufftComplex *)d_data_spec) == CUFFT_SUCCESS);
    modulateAndNormalize3D(d_data_spec, d_kernel_spec, fftZ, fftY, fftX, /*padding=*/1);
    cudaDeviceSynchronize();
    REQUIRE(cufftExecC2R(planC2R, (cufftComplex *)d_data_spec, (cufftReal *)d_data_pad) == CUFFT_SUCCESS);
    cudaDeviceSynchronize();

    unpadResult3D(d_out, d_data_pad, dZ, dY, dX, fftZ, fftY, fftX);
    cudaDeviceSynchronize();

    cufftDestroy(planR2C);
    cufftDestroy(planC2R);

    std::vector<float> h_out(h_data.size(), 0.0f);
    cudaMemcpy(h_out.data(), d_out, sizeof(float) * h_out.size(), cudaMemcpyDeviceToHost);

    cudaFree(d_data);
    cudaFree(d_kernel);
    cudaFree(d_data_pad);
    cudaFree(d_kernel_pad);
    cudaFree(d_data_spec);
    cudaFree(d_kernel_spec);
    cudaFree(d_out);

    // Compare only interior voxels: [1..n-2] in each dimension.
    double num = 0.0, den = 0.0;
    double max_abs = 0.0;
    for (int z = 1; z < nz - 1; z++)
        for (int y = 1; y < ny - 1; y++)
            for (int x = 1; x < nx - 1; x++)
            {
                const int i = x + nx * (y + ny * z);
                const double ref = cpu_out[i];
                const double got = (double)h_out[i];
                const double diff = got - ref;
                num += diff * diff;
                den += ref * ref;
                max_abs = std::max(max_abs, std::fabs(diff));
            }
    const double rel_l2 = (den > 0.0) ? std::sqrt(num / den) : std::sqrt(num);
    INFO("rel_l2=" << rel_l2 << " max_abs=" << max_abs);
    CHECK(rel_l2 < 1e-4);
    CHECK(max_abs < 1e-3);
#endif
}
