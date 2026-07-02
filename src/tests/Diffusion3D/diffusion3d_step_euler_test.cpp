/**
 * @file diffusion3d_step_euler_test.cpp
 * @brief Catch2 tests for diffusion3d_step_euler_scalar / diffusion3d_step_euler_gpu.
 */

#include <cmath>

#include <catch2/catch.hpp>

#include "diffusion3d_step_euler_cpu.h"
#include "diffusion3d_gpu.h"
#include "scalar_field_grid.h"

#ifdef DIFFUSION3D_CUDA
#include <cuda_runtime.h>
#endif

namespace
{

bool nearly_equal_vec(const std::vector<double> &a, const std::vector<double> &b, double eps)
{
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        if (std::fabs(a[i] - b[i]) > eps)
            return false;
    }
    return true;
}

void fill_sin_ramp(ScalarFieldGrid &grid)
{
    for (int z = 0; z < grid.nz_; ++z)
    {
        for (int y = 0; y < grid.ny_; ++y)
        {
            for (int x = 0; x < grid.nx_; ++x)
            {
                grid.at(x, y, z) =
                    std::sin(0.3 * static_cast<double>(x)) + 0.2 * static_cast<double>(y + z);
            }
        }
    }
}

void cpu_advance(ScalarFieldGrid &grid, double h, double D, int n_steps, double dt)
{
    for (int i = 0; i < n_steps; ++i)
        diffusion3d_step_euler_scalar(grid, h, D, dt);
}

} // namespace

TEST_CASE("diffusion3d_step_euler_scalar identical ICs stay aligned over multiple steps", "[cpu][euler]")
{
    const int nx = 6, ny = 5, nz = 4;
    const double h = 0.4, D = 0.08, dt = 0.02;
    const int nsteps = 4;

    ScalarFieldGrid a(nx, ny, nz);
    ScalarFieldGrid b(nx, ny, nz);
    fill_sin_ramp(a);
    fill_sin_ramp(b);

    cpu_advance(a, h, D, nsteps, dt);
    cpu_advance(b, h, D, nsteps, dt);

    CHECK(nearly_equal_vec(a.data_, b.data_, 1e-13));
}

#ifdef DIFFUSION3D_CUDA

TEST_CASE("diffusion3d_step_euler_gpu chained matches CPU euler chain", "[gpu][euler]")
{
    const int nx = 6, ny = 5, nz = 4;
    const double h = 0.4, D = 0.08, dt = 0.02;
    const int nsteps = 4;

    ScalarFieldGrid cpu_grid(nx, ny, nz);
    ScalarFieldGrid gpu_grid(nx, ny, nz);
    fill_sin_ramp(cpu_grid);
    fill_sin_ramp(gpu_grid);

    cpu_advance(cpu_grid, h, D, nsteps, dt);

    double *d_u = nullptr;
    double *d_u_next = nullptr;
    const size_t bytes = static_cast<size_t>(gpu_grid.size()) * sizeof(double);
    REQUIRE(cudaMalloc(&d_u, bytes) == cudaSuccess);
    REQUIRE(cudaMalloc(&d_u_next, bytes) == cudaSuccess);
    REQUIRE(cudaMemcpy(d_u, gpu_grid.data_.data(), bytes, cudaMemcpyHostToDevice) == cudaSuccess);

    for (int s = 0; s < nsteps; ++s)
        diffusion3d_step_euler_gpu(nx, ny, nz, h, D, dt, d_u, d_u_next);

    REQUIRE(cudaMemcpy(gpu_grid.data_.data(), d_u, bytes, cudaMemcpyDeviceToHost) == cudaSuccess);
    cudaFree(d_u);
    cudaFree(d_u_next);

    CHECK(nearly_equal_vec(cpu_grid.data_, gpu_grid.data_, 1e-13));
}

#endif // DIFFUSION3D_CUDA
