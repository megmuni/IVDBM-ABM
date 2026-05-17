/**
 * @file diffusion3d_solver.cpp
 * @brief Implements DiffusionSolver (substepping + CPU/GPU stencil / GPU FFT dispatch).
 */

#include "diffusion3d_solver.h"
#include "diffusion3d_timestep.h"

#ifdef GPU_DIFFUSE
#include <cuda_runtime.h>
#include "diffusion3d_cuda_check.h"
#include "diffusion3d_fft_scratch.h"
#endif

namespace
{

DiffusionBackend resolve_backend(DiffusionBackend b)
{
#ifdef GPU_DIFFUSE
    // Prefer the FFT path when CUDA is available: one cufft round-trip per tick is cheaper
    // than n_sub explicit Euler launches once n_sub grows.
    if (b == DiffusionBackend::Auto)
        return DiffusionBackend::GpuFftPrecomputed;
#else
    if (b == DiffusionBackend::Auto)
        return DiffusionBackend::CpuStencil;
#endif
    return b;
}

#ifdef GPU_DIFFUSE
/** Smallest multiple of `align` that is >= `a` (align > 0). */
inline int align_dimension_up(int a, int align)
{
    assert(align > 0);
    const int r = a % align;
    return (r != 0) ? (a - r + align) : a;
}

/** Even R2C real-axis length: use user_real_extent if nonzero, else derived padding (min 32, step 16). */
int derive_fft_axis_len(int domain_len, int user_real_extent)
{
    assert(domain_len > 0);
    if (user_real_extent > 0)
    {
        assert(user_real_extent % 2 == 0);
        return user_real_extent;
    }

    const int need = domain_len + 2;
    int v = align_dimension_up(need, 16);
    if (v < 32)
        v = 32;
    if (v % 2 != 0)
        ++v;
    return v;
}
#endif

} // namespace

DiffusionSolver::DiffusionSolver(DiffusionParams params)
    : params_(params),
      resolved_backend_(resolve_backend(params.backend))
{
    assert(params_.safety > 0.0 && params_.safety <= 1.0);
#ifndef GPU_DIFFUSE
    assert(resolved_backend_ != DiffusionBackend::GpuStencil &&
           "DiffusionBackend::GpuStencil requires -DGPU_DIFFUSE=ON");
    assert(resolved_backend_ != DiffusionBackend::GpuFftPrecomputed &&
           "DiffusionBackend::GpuFftPrecomputed requires -DGPU_DIFFUSE=ON");
#endif
}

DiffusionSolver::~DiffusionSolver()
{
    free_gpu_buffers();
#ifdef GPU_DIFFUSE
    diffusion3d_fft_scratch_destroy(fft_scratch_);
    fft_scratch_ = nullptr;
#endif
}

void DiffusionSolver::ensure_fft_scratch(Diffusion3DContext &ctx)
{
#ifdef GPU_DIFFUSE
    const int fx = derive_fft_axis_len(ctx.nx, params_.fft_real_extent_x);
    const int fy = derive_fft_axis_len(ctx.ny, params_.fft_real_extent_y);
    const int fz = derive_fft_axis_len(ctx.nz, params_.fft_real_extent_z);

    if (fft_scratch_ != nullptr && fft_cached_nx_ == ctx.nx && fft_cached_ny_ == ctx.ny &&
        fft_cached_nz_ == ctx.nz && fft_cached_fx_ == fx && fft_cached_fy_ == fy &&
        fft_cached_fz_ == fz)
        return;

    diffusion3d_fft_scratch_destroy(fft_scratch_);
    fft_scratch_ =
        diffusion3d_fft_scratch_create(ctx.nx, ctx.ny, ctx.nz, fx, fy, fz);
    fft_cached_nx_ = ctx.nx;
    fft_cached_ny_ = ctx.ny;
    fft_cached_nz_ = ctx.nz;
    fft_cached_fx_ = fx;
    fft_cached_fy_ = fy;
    fft_cached_fz_ = fz;
#else
    (void)ctx;
#endif
}

void DiffusionSolver::ensure_gpu_buffers(std::size_t n)
{
#ifdef GPU_DIFFUSE
    if (n == 0)
        return;
    if (gpu_n_ == n && d_u_ != nullptr && d_u_next_ != nullptr)
        return;

    free_gpu_buffers();

    const std::size_t bytes = n * sizeof(double);
    cuda_ok(cudaMalloc(reinterpret_cast<void **>(&d_u_), bytes), "cudaMalloc d_u_");
    cuda_ok(cudaMalloc(reinterpret_cast<void **>(&d_u_next_), bytes), "cudaMalloc d_u_next_");
    gpu_n_ = n;
#else
    (void)n;
#endif
}

void DiffusionSolver::free_gpu_buffers()
{
#ifdef GPU_DIFFUSE
    if (d_u_ != nullptr)
    {
        cuda_ok(cudaFree(d_u_), "cudaFree d_u_");
        d_u_ = nullptr;
    }
    if (d_u_next_ != nullptr)
    {
        cuda_ok(cudaFree(d_u_next_), "cudaFree d_u_next_");
        d_u_next_ = nullptr;
    }
    gpu_n_ = 0;
#endif
}

void DiffusionSolver::configure_tick(Diffusion3DContext &ctx, double tick_dt)
{
    assert(tick_dt >= 0.0);

    cfg_nx_ = ctx.nx;
    cfg_ny_ = ctx.ny;
    cfg_nz_ = ctx.nz;
    cfg_h_ = ctx.h;
    cfg_D_ = ctx.D;

    tick_dt_ = tick_dt;

    const double dt_max = compute_stability_constraint(ctx.h, ctx.D, params_.safety);
    const SubstepPlan plan = plan_substeps(tick_dt_, dt_max);
    n_sub_ = plan.n_sub;
    dt_sub_ = plan.dt_sub;

    if (resolved_backend_ == DiffusionBackend::GpuStencil)
    {
#ifdef GPU_DIFFUSE
        ensure_gpu_buffers(static_cast<std::size_t>(ctx.n));
#endif
    }

    if (resolved_backend_ == DiffusionBackend::GpuFftPrecomputed)
    {
#ifdef GPU_DIFFUSE
        ensure_fft_scratch(ctx);
        // Precompute (K_hat_one_step)^n_sub for this fixed tick.
        diffusion3d_fft_scratch_rebuild_operator(fft_scratch_, ctx.h, ctx.D, tick_dt_, params_.safety);
#endif
    }

    configured_ = true;
}

void DiffusionSolver::advance_tick(Diffusion3DContext &ctx)
{
    assert(configured_);
    assert(ctx.nx == cfg_nx_ && ctx.ny == cfg_ny_ && ctx.nz == cfg_nz_);
    assert(ctx.h == cfg_h_);
    assert(ctx.D == cfg_D_);

    if (resolved_backend_ == DiffusionBackend::GpuFftPrecomputed)
    {
#ifdef GPU_DIFFUSE
        ensure_fft_scratch(ctx);
        diffusion3d_fft_scratch_apply(fft_scratch_, ctx);
#endif
        return;
    }

    if (resolved_backend_ == DiffusionBackend::GpuStencil)
    {
#ifdef GPU_DIFFUSE
        ensure_gpu_buffers(static_cast<std::size_t>(ctx.n));
        assert(static_cast<std::size_t>(ctx.n) == gpu_n_);
        const std::size_t bytes = gpu_n_ * sizeof(double);
        cuda_ok(cudaMemcpy(d_u_, ctx.u.data(), bytes, cudaMemcpyHostToDevice), "H2D ctx.u");

        double *d_u = d_u_;
        double *d_next = d_u_next_;
        for (int i = 0; i < n_sub_; ++i)
            diffusion3d_step_euler_gpu(ctx, dt_sub_, d_u, d_next);

        d_u_ = d_u;
        d_u_next_ = d_next;

        cuda_ok(cudaMemcpy(ctx.u.data(), d_u_, bytes, cudaMemcpyDeviceToHost), "D2H ctx.u");
#else
        assert(false && "GpuStencil requires GPU_DIFFUSE build");
#endif
        return;
    }

    for (int i = 0; i < n_sub_; ++i)
        diffusion3d_step_euler_cpu(ctx, dt_sub_);
}
