/**
 * @file diffusion3d_step_euler_test.cpp
 * @brief Catch2 tests for diffusion3d_step_euler_cpu / diffusion3d_step_euler_gpu.
 */

#include <cmath>

#include <catch2/catch.hpp>

#include "diffusion3d_common.h"

#ifdef GPU_DIFFUSE
#include <cuda_runtime.h>
#endif

static bool nearly_equal_vec(const std::vector<double> &a, const std::vector<double> &b,
                             double eps)
{
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        if (std::fabs(a[i] - b[i]) > eps)
        {
            return false;
        }
    }
    return true;
}

/** @brief Advance `ctx` by `n_steps` explicit Euler steps at fixed `dt`. */
static void cpu_advance(Diffusion3DContext &ctx, int n_steps, double dt)
{
    for (int i = 0; i < n_steps; ++i)
        diffusion3d_step_euler_cpu(ctx, dt);
}

TEST_CASE("diffusion3d_step_euler_cpu identical ICs stay aligned over multiple steps", "[cpu][euler]")
{
    const int nx = 6, ny = 5, nz = 4;
    const double h = 0.4, D = 0.08, dt = 0.02;
    const int nsteps = 4;

    Diffusion3DContext a = Diffusion3DContext::make(nx, ny, nz, h, D);
    Diffusion3DContext b = Diffusion3DContext::make(nx, ny, nz, h, D);

    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const double v =
                    std::sin(0.3 * static_cast<double>(x)) + 0.2 * static_cast<double>(y + z);
                a.set_at_coord(x, y, z, v);
                b.set_at_coord(x, y, z, v);
            }
        }
    }

    cpu_advance(a, nsteps, dt);
    cpu_advance(b, nsteps, dt);

    CHECK(nearly_equal_vec(a.u, b.u, 1e-13));
}

#ifdef GPU_DIFFUSE

TEST_CASE("diffusion3d_step_euler_gpu chained matches CPU euler chain", "[gpu][euler]")
{
    const int nx = 6, ny = 5, nz = 4;
    const double h = 0.4, D = 0.08, dt = 0.02;
    const int nsteps = 4;

    Diffusion3DContext cpu_ctx = Diffusion3DContext::make(nx, ny, nz, h, D);
    Diffusion3DContext gpu_ctx = Diffusion3DContext::make(nx, ny, nz, h, D);

    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const double v =
                    std::sin(0.3 * static_cast<double>(x)) + 0.2 * static_cast<double>(y + z);
                cpu_ctx.set_at_coord(x, y, z, v);
                gpu_ctx.set_at_coord(x, y, z, v);
            }
        }
    }

    cpu_advance(cpu_ctx, nsteps, dt);

    double *d_u = nullptr;
    double *d_u_next = nullptr;
    const size_t bytes = static_cast<size_t>(gpu_ctx.n) * sizeof(double);
    REQUIRE(cudaMalloc(&d_u, bytes) == cudaSuccess);
    REQUIRE(cudaMalloc(&d_u_next, bytes) == cudaSuccess);
    REQUIRE(cudaMemcpy(d_u, gpu_ctx.u.data(), bytes, cudaMemcpyHostToDevice) == cudaSuccess);

    for (int s = 0; s < nsteps; ++s)
        diffusion3d_step_euler_gpu(gpu_ctx, dt, d_u, d_u_next);

    REQUIRE(cudaMemcpy(gpu_ctx.u.data(), d_u, bytes, cudaMemcpyDeviceToHost) == cudaSuccess);
    cudaFree(d_u);
    cudaFree(d_u_next);

    CHECK(nearly_equal_vec(cpu_ctx.u, gpu_ctx.u, 1e-13));
}

#endif // GPU_DIFFUSE
