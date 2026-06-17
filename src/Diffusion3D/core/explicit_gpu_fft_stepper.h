#ifndef DIFFUSION3D_EXPLICIT_GPU_FFT_STEPPER_H
#define DIFFUSION3D_EXPLICIT_GPU_FFT_STEPPER_H

/**
 * @file explicit_gpu_fft_stepper.h
 * @brief GPU FFT precomputed multi-species stepper.
 */

#include <map>
#include <vector>
#include "multi_species_diffusion_engine.h"
#include "diffusion3d_timestep.h"

struct DiffusionFftScratch;

class ExplicitGpuFftStepper : public MultiSpeciesDiffusionEngine
{
public:
    ExplicitGpuFftStepper();
    ~ExplicitGpuFftStepper() override;

    ExplicitGpuFftStepper(const ExplicitGpuFftStepper &) = delete;
    ExplicitGpuFftStepper &operator=(const ExplicitGpuFftStepper &) = delete;

    void configure_species_interval(
        const MultiSpeciesFieldGrid &grid,
        const MultiSpeciesDiffusionSettings &settings,
        double tick_dt) override;

    void advance_species_interval(MultiSpeciesFieldGrid &grid) override;

    std::vector<SpeciesId> supported_species() const override;

    DiffusionAlgorithm resolved_algorithm() const override
    {
        return DiffusionAlgorithm::GpuFftPrecomputed;
    }

private:
    void destroy_fft_scratch(SpeciesId id);

    bool configured_ = false;
    int cfg_nx_ = 0, cfg_ny_ = 0, cfg_nz_ = 0;
    double cfg_h_ = 1.0;

    MultiSpeciesDiffusionSettings settings_;
    std::vector<SpeciesId> species_ids_;
    std::map<SpeciesId, DiffusionFftScratch *> fft_scratch_;
};

#endif // DIFFUSION3D_EXPLICIT_GPU_FFT_STEPPER_H
