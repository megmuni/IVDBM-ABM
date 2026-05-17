/**
 * @file diffusion3d_solver_test.cpp
 * @brief Catch2 tests for DiffusionSolver (substeps + backends).
 */

#include <cmath>

#include <catch2/catch.hpp>

#include "diffusion3d_common.h"
#include "diffusion3d_solver.h"
#include "diffusion3d_timestep.h"

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

} // namespace

TEST_CASE("DiffusionSolver CpuStencil splits tick into stable substeps", "[solver][cpu]")
{
    const int nx = 5, ny = 4, nz = 3;
    const double h = 1.0;
    const double D = 1.0;
    const double safety = 1.0;
    const double tick_dt = 0.5;

    const double dt_max = compute_stability_constraint(h, D, safety);
    const SubstepPlan plan = plan_substeps(tick_dt, dt_max);
    REQUIRE(plan.n_sub >= 2);

    DiffusionParams params;
    params.safety = safety;
    params.backend = DiffusionBackend::CpuStencil;
    DiffusionSolver solver(params);
    REQUIRE(solver.resolved_backend() == DiffusionBackend::CpuStencil);

    Diffusion3DContext ref = Diffusion3DContext::make(nx, ny, nz, h, D);
    Diffusion3DContext sol = Diffusion3DContext::make(nx, ny, nz, h, D);
    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const double v = 0.01 * static_cast<double>(x + y * nx + z * nx * ny);
                ref.set_at_coord(x, y, z, v);
                sol.set_at_coord(x, y, z, v);
            }
        }
    }

    for (int i = 0; i < plan.n_sub; ++i)
        diffusion3d_step_euler_cpu(ref, plan.dt_sub);

    solver.configure_tick(sol, tick_dt);
    solver.advance_tick(sol);

    CHECK(nearly_equal_vec(ref.u, sol.u, 1e-13));
}

TEST_CASE("DiffusionSolver Auto resolves per build (CpuStencil w/o CUDA, GpuFftPrecomputed w/ CUDA)",
          "[solver][backend]")
{
#ifndef GPU_DIFFUSE
    DiffusionSolver s(DiffusionParams{});
    CHECK(s.resolved_backend() == DiffusionBackend::CpuStencil);
#else
    DiffusionSolver s(DiffusionParams{});
    CHECK(s.resolved_backend() == DiffusionBackend::GpuFftPrecomputed);
#endif
}

#ifdef GPU_DIFFUSE

TEST_CASE("DiffusionSolver GpuStencil matches CpuStencil micro-steps", "[solver][gpu]")
{
    const int nx = 6, ny = 5, nz = 4;
    const double h = 0.4, D = 0.08;
    const double safety = 1.0;
    const double tick_dt = 0.07;

    const double dt_max = compute_stability_constraint(h, D, safety);
    const SubstepPlan plan = plan_substeps(tick_dt, dt_max);
    REQUIRE(plan.n_sub >= 1);

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

    for (int i = 0; i < plan.n_sub; ++i)
        diffusion3d_step_euler_cpu(cpu_ctx, plan.dt_sub);

    DiffusionParams p;
    p.safety = safety;
    p.backend = DiffusionBackend::GpuStencil;
    DiffusionSolver gpu_solver(p);
    gpu_solver.configure_tick(gpu_ctx, tick_dt);
    gpu_solver.advance_tick(gpu_ctx);

    CHECK(nearly_equal_vec(cpu_ctx.u, gpu_ctx.u, 1e-12));
}

TEST_CASE("DiffusionSolver GpuFftPrecomputed matches CpuStencil in interior (small tick)", "[solver][gpu][fft]")
{
    const int nx = 6, ny = 6, nz = 6;
    const double tick_dt = 1e-3;
    const double h = 1.0;
    const double D = 0.1;
    const double safety = 1.0;

    Diffusion3DContext cpu_ctx = Diffusion3DContext::make(nx, ny, nz, h, D);
    Diffusion3DContext fft_ctx = Diffusion3DContext::make(nx, ny, nz, h, D);

    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const double v =
                    0.01 * static_cast<double>((z + 1) * 100 + (y + 1) * 10 + (x + 1));
                cpu_ctx.set_at_coord(x, y, z, v);
                fft_ctx.set_at_coord(x, y, z, v);
            }
        }
    }

    DiffusionParams cpu_p;
    cpu_p.backend = DiffusionBackend::CpuStencil;
    cpu_p.safety = safety;
    DiffusionSolver cpu_solver(cpu_p);
    cpu_solver.configure_tick(cpu_ctx, tick_dt);
    cpu_solver.advance_tick(cpu_ctx);

    DiffusionParams fft_p;
    fft_p.backend = DiffusionBackend::GpuFftPrecomputed;
    fft_p.safety = safety;
    fft_p.fft_real_extent_x = 8;
    fft_p.fft_real_extent_y = 8;
    fft_p.fft_real_extent_z = 8;
    DiffusionSolver fft_solver(fft_p);
    fft_solver.configure_tick(fft_ctx, tick_dt);
    fft_solver.advance_tick(fft_ctx);

    double max_abs = 0.0;
    double num = 0.0;
    double den = 0.0;
    for (int z = 1; z < nz - 1; ++z)
    {
        for (int y = 1; y < ny - 1; ++y)
        {
            for (int x = 1; x < nx - 1; ++x)
            {
                const int i = cpu_ctx.idx(x, y, z);
                const double diff = fft_ctx.u[static_cast<size_t>(i)] - cpu_ctx.u[static_cast<size_t>(i)];
                max_abs = std::max(max_abs, std::fabs(diff));
                num += diff * diff;
                den += cpu_ctx.u[static_cast<size_t>(i)] * cpu_ctx.u[static_cast<size_t>(i)];
            }
        }
    }
    const double rel_l2 = (den > 0.0) ? std::sqrt(num / den) : std::sqrt(num);
    INFO("rel_l2=" << rel_l2 << " max_abs=" << max_abs);
    CHECK(rel_l2 < 1e-4);
    CHECK(max_abs < 1e-3);
}

#endif // GPU_DIFFUSE
