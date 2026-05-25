/**
 * @file make_multi_species_diffusion_engine.h
 * @brief Factory for creating multi-species diffusion engines.
 */

#ifndef DIFFUSION3D_MAKE_MULTI_SPECIES_DIFFUSION_ENGINE_H
#define DIFFUSION3D_MAKE_MULTI_SPECIES_DIFFUSION_ENGINE_H

#include <memory>
#include <stdexcept>
#include "multi_species_diffusion_engine.h"
#include "explicit_multi_species_heat_stepper.h"
#include "explicit_gpu_stencil_stepper.h"
#include "explicit_gpu_fft_stepper.h"

inline DiffusionAlgorithm ResolveDiffusionAlgorithm(const MultiSpeciesDiffusionSettings &settings)
{
#ifndef DIFFUSION3D_CUDA
    (void)settings;
    return DiffusionAlgorithm::ExplicitHeatEquation;
#else
    return settings.algorithm;
#endif
}

/**
 * @brief Create a MultiSpeciesDiffusionEngine for the requested algorithm.
 *
 * Caller must invoke configure_species_interval(grid, settings, tick_dt) before advance.
 */
inline std::unique_ptr<MultiSpeciesDiffusionEngine> MakeMultiSpeciesDiffusionEngine(
    const MultiSpeciesDiffusionSettings &settings)
{
    settings.validate();

    switch (ResolveDiffusionAlgorithm(settings))
    {
    case DiffusionAlgorithm::ExplicitHeatEquation:
        return std::make_unique<ExplicitMultiSpeciesHeatStepper>();
    case DiffusionAlgorithm::GpuStencil:
        return std::make_unique<ExplicitGpuStencilStepper>();
    case DiffusionAlgorithm::GpuFftPrecomputed:
        return std::make_unique<ExplicitGpuFftStepper>();
    default:
        break;
    }
    throw std::logic_error("unsupported DiffusionAlgorithm");
}

#endif // DIFFUSION3D_MAKE_MULTI_SPECIES_DIFFUSION_ENGINE_H
