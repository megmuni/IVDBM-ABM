#ifndef DIFFUSION3D_EXPLICIT_GPU_STENCIL_STEPPER_H
#define DIFFUSION3D_EXPLICIT_GPU_STENCIL_STEPPER_H

/**
 * @file explicit_gpu_stencil_stepper.h
 * @brief GPU 6-point stencil multi-species stepper.
 */

#include <map>
#include <vector>
#include "multi_species_diffusion_engine.h"
#include "diffusion3d_timestep.h"

class ExplicitGpuStencilStepper : public MultiSpeciesDiffusionEngine
{
public:
    ExplicitGpuStencilStepper();
    ~ExplicitGpuStencilStepper() override;

    ExplicitGpuStencilStepper(const ExplicitGpuStencilStepper &) = delete;
    ExplicitGpuStencilStepper &operator=(const ExplicitGpuStencilStepper &) = delete;

    void configure_species_interval(
        const MultiSpeciesFieldGrid &grid,
        const MultiSpeciesDiffusionSettings &settings,
        double tick_dt) override;

    void advance_species_interval(MultiSpeciesFieldGrid &grid) override;

    std::vector<SpeciesId> supported_species() const override;

    DiffusionAlgorithm resolved_algorithm() const override
    {
        return DiffusionAlgorithm::GpuStencil;
    }

private:
    void free_gpu_buffers();
    void ensure_gpu_buffers(std::size_t n);

    bool configured_ = false;
    double tick_dt_ = 0.0;
    int cfg_nx_ = 0, cfg_ny_ = 0, cfg_nz_ = 0;
    double cfg_h_ = 1.0;

    MultiSpeciesDiffusionSettings settings_;
    std::vector<SpeciesId> species_ids_;
    std::map<SpeciesId, SubstepPlan> plan_per_species_;

    double *d_u_ = nullptr;
    double *d_u_next_ = nullptr;
    std::size_t gpu_n_ = 0;
};

#endif // DIFFUSION3D_EXPLICIT_GPU_STENCIL_STEPPER_H
