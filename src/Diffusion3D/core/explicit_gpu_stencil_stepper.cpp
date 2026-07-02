/**
 * @file explicit_gpu_stencil_stepper.cpp
 * @brief GPU stencil multi-species stepper.
 */

#include "explicit_gpu_stencil_stepper.h"
#include "diffusion3d_gpu.h"
#include <stdexcept>

#ifdef DIFFUSION3D_CUDA
#include <cuda_runtime.h>
#include "diffusion3d_cuda_check.h"
#endif

ExplicitGpuStencilStepper::ExplicitGpuStencilStepper() = default;

ExplicitGpuStencilStepper::~ExplicitGpuStencilStepper()
{
    free_gpu_buffers();
}

void ExplicitGpuStencilStepper::free_gpu_buffers()
{
#ifdef DIFFUSION3D_CUDA
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

void ExplicitGpuStencilStepper::ensure_gpu_buffers(std::size_t n)
{
#ifdef DIFFUSION3D_CUDA
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

void ExplicitGpuStencilStepper::configure_species_interval(
    const MultiSpeciesFieldGrid &grid,
    const MultiSpeciesDiffusionSettings &settings,
    double tick_dt)
{
#ifndef DIFFUSION3D_CUDA
    throw std::runtime_error("ExplicitGpuStencilStepper requires DIFFUSION3D_CUDA=ON");
#else
    settings.validate();
    species_ids_ = grid.species();
    if (species_ids_.empty())
        throw std::runtime_error("Grid must contain at least one species");

    settings_ = settings;
    cfg_nx_ = grid.nx_;
    cfg_ny_ = grid.ny_;
    cfg_nz_ = grid.nz_;
    cfg_h_ = grid.h_;
    tick_dt_ = tick_dt;

    plan_per_species_.clear();
    for (auto id : species_ids_)
    {
        const double D = settings.species_diffusivities.at(id);
        const double dt_max = compute_stability_constraint(cfg_h_, D, settings.safety);
        plan_per_species_[id] = plan_substeps(tick_dt_, dt_max);
    }

    ensure_gpu_buffers(static_cast<std::size_t>(cfg_nx_) * cfg_ny_ * cfg_nz_);
    configured_ = true;
#endif
}

void ExplicitGpuStencilStepper::advance_species_interval(MultiSpeciesFieldGrid &grid)
{
#ifndef DIFFUSION3D_CUDA
    throw std::runtime_error("ExplicitGpuStencilStepper requires DIFFUSION3D_CUDA=ON");
#else
    if (!configured_)
        throw std::runtime_error("Stepper not configured; call configure_species_interval() first");

    if (grid.nx_ != cfg_nx_ || grid.ny_ != cfg_ny_ || grid.nz_ != cfg_nz_ || grid.h_ != cfg_h_)
        throw std::runtime_error("Grid dimensions changed since configure_species_interval()");

    const std::size_t bytes = gpu_n_ * sizeof(double);

    for (auto species_id : species_ids_)
    {
        ScalarFieldGrid &species_grid = grid.grid(species_id);
        const double D = settings_.species_diffusivities.at(species_id);
        const SubstepPlan &plan = plan_per_species_.at(species_id);

        cuda_ok(cudaMemcpy(d_u_, species_grid.data_.data(), bytes, cudaMemcpyHostToDevice), "H2D u");

        double *d_u = d_u_;
        double *d_u_next = d_u_next_;
        for (int i = 0; i < plan.n_sub; ++i)
            diffusion3d_step_euler_gpu(cfg_nx_, cfg_ny_, cfg_nz_, cfg_h_, D, plan.dt_sub, d_u, d_u_next);

        d_u_ = d_u;
        d_u_next_ = d_u_next;

        cuda_ok(cudaMemcpy(species_grid.data_.data(), d_u_, bytes, cudaMemcpyDeviceToHost), "D2H u");
    }
#endif
}

std::vector<SpeciesId> ExplicitGpuStencilStepper::supported_species() const
{
    return species_ids_;
}
