#include <cassert>
#include <filesystem>
#include <system_error>
#include <iostream>
#include <string>

#include "diffusion3d_common.h"
#include "diffusion3d_solver.h"
#include "logger.h"
#include "vtk_export.h"

/**
 * @file main.cpp
 * @brief Entrypoint: CPU/GPU stencil + optional GPU FFT (precomputed (K_hat)^N per advance) writing ParaView `.vti`.
 */

int main()
{
    const std::string cpu_outdir = "cpu_series";
    const std::string gpu_outdir = "gpu_series";
    const std::string fft_outdir = "fft_series";

    // Ensure ParaView output dirs exist (relative to current working directory).
    // Use error_code: create_directories returns false if the path already exists, which is OK.
    std::error_code dir_ec;
    std::filesystem::create_directories(cpu_outdir, dir_ec);
    assert(!dir_ec);
    std::filesystem::create_directories(gpu_outdir, dir_ec);
    assert(!dir_ec);
    std::filesystem::create_directories(fft_outdir, dir_ec);
    assert(!dir_ec);

    // Test scenario:
    // - user picks macro tick (`tick_dt`), grid spacing (`h`), and diffusivity (`D`)
    // - solver enforces explicit-Euler stability internally by choosing (n_sub, dt_sub)
    const int nsteps = 4;
    const int nx = 8, ny = 8, nz = 8;
    const double tick_dt = 2.0;
    const double h = 1.0;
    const double D = 0.1;

    // --- CPU TIME SERIES ---
    STEP_LOG("main: creating CPU context");
    Diffusion3DContext cpu_ctx = Diffusion3DContext::make_with_reference_dt(nx, ny, nz, h, D, tick_dt);
    cpu_ctx.set_at_coord(nx / 2, ny / 2, nz / 2, 1.0);

    DiffusionParams cpu_params;
    cpu_params.backend = DiffusionBackend::CpuStencil;
    DiffusionSolver cpu_solver(cpu_params);
    cpu_solver.configure_tick(cpu_ctx, tick_dt);
    std::cout << "main: CPU configured tick_dt=" << cpu_solver.tick_dt()
              << " n_sub=" << cpu_solver.n_sub()
              << " dt_sub=" << cpu_solver.dt_sub() << "\n";
    assert(cpu_solver.n_sub() > 1 && "pick tick_dt/h/D that triggers solver substepping");

    // Export t = 0 … nsteps so file `t{k}` is the state after k steps (aligned with fft_series naming).
    for (int t = 0; t <= nsteps; ++t)
    {
        STEP_LOG("main: CPU snapshot t=" + std::to_string(t));
        char fname[256];
        snprintf(fname, sizeof(fname), "%s/diffusion3d_t%d.vti", cpu_outdir.c_str(), t);
        export_field_to_vti(cpu_ctx.u, nx, ny, nz, fname, "u");
        if (t < nsteps)
            cpu_solver.advance_tick(cpu_ctx);
    }

    // --- GPU STENCIL TIME SERIES ---
#ifdef GPU_DIFFUSE
    STEP_LOG("main: creating GPU stencil context");
    Diffusion3DContext gpu_ctx = Diffusion3DContext::make_with_reference_dt(nx, ny, nz, h, D, tick_dt);
    gpu_ctx.set_at_coord(nx / 2, ny / 2, nz / 2, 1.0);

    DiffusionParams gpu_params;
    gpu_params.backend = DiffusionBackend::GpuStencil;
    DiffusionSolver gpu_solver(gpu_params);
    gpu_solver.configure_tick(gpu_ctx, tick_dt);
    std::cout << "main: GPU stencil configured tick_dt=" << gpu_solver.tick_dt()
              << " n_sub=" << gpu_solver.n_sub()
              << " dt_sub=" << gpu_solver.dt_sub() << "\n";
    assert(gpu_solver.n_sub() > 1 && "pick tick_dt/h/D that triggers solver substepping");

    for (int t = 0; t <= nsteps; ++t)
    {
        STEP_LOG("main: GPU stencil snapshot t=" + std::to_string(t));
        char fname[256];
        snprintf(fname, sizeof(fname), "%s/diffusion3d_t%d.vti", gpu_outdir.c_str(), t);
        export_field_to_vti(gpu_ctx.u, nx, ny, nz, fname, "u");
        if (t < nsteps)
            gpu_solver.advance_tick(gpu_ctx);
    }
#else
    STEP_LOG("main: GPU stencil series — CUDA off; exporting initial field each frame");
    Diffusion3DContext gpu_ctx = Diffusion3DContext::make_with_reference_dt(nx, ny, nz, h, D, tick_dt);
    gpu_ctx.set_at_coord(nx / 2, ny / 2, nz / 2, 1.0);

    for (int t = 0; t <= nsteps; ++t)
    {
        STEP_LOG("main: GPU stencil snapshot t=" + std::to_string(t));
        char fname[256];
        snprintf(fname, sizeof(fname), "%s/diffusion3d_t%d.vti", gpu_outdir.c_str(), t);
        export_field_to_vti(gpu_ctx.u, nx, ny, nz, fname, "u");
    }
#endif

#ifdef GPU_DIFFUSE
    // --- GPU FFT: one padded conv per tick via diffusion3d_fft_scratch (composed K_hat^n_sub) ---
    STEP_LOG("main: FFT path — GpuFftPrecomputed");
    Diffusion3DContext fft_ctx = Diffusion3DContext::make_with_reference_dt(nx, ny, nz, h, D, tick_dt);
    fft_ctx.set_at_coord(nx / 2, ny / 2, nz / 2, 1.0);

    DiffusionParams fft_params;
    fft_params.backend = DiffusionBackend::GpuFftPrecomputed;
    fft_params.fft_real_extent_x = 32;
    fft_params.fft_real_extent_y = 32;
    fft_params.fft_real_extent_z = 32;
    DiffusionSolver fft_solver(fft_params);
    fft_solver.configure_tick(fft_ctx, tick_dt);
    std::cout << "main: FFT configured tick_dt=" << fft_solver.tick_dt()
              << " n_sub=" << fft_solver.n_sub()
              << " dt_sub=" << fft_solver.dt_sub() << "\n";
    assert(fft_solver.n_sub() > 1 && "pick tick_dt/h/D that triggers solver substepping");

    {
        char fname[256];
        snprintf(fname, sizeof(fname), "%s/diffusion3d_t0.vti", fft_outdir.c_str());
        export_to_vti(fft_ctx, fname);
    }

    for (int iter = 0; iter < nsteps; ++iter)
        fft_solver.advance_tick(fft_ctx);

    {
        char fname[256];
        snprintf(fname, sizeof(fname), "%s/diffusion3d_t%d.vti", fft_outdir.c_str(), nsteps);
        export_to_vti(fft_ctx, fname);
    }
    std::cout << "main: FFT series completed" << std::endl;
#else
    std::cout << "main: GPU_DIFFUSE=OFF — skip FFT convolution series (build with -DGPU_DIFFUSE=ON for fft_series/)\n";
#endif

    return 0;
}
