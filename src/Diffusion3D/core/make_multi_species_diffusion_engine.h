/**
 * @file make_multi_species_diffusion_engine.h
 * @brief Factory for creating multi-species diffusion engines for Diffusion3D.
 *
 * Provides MakeMultiSpeciesDiffusionEngine for constructing engine instances.
 */

#ifndef DIFFUSION3D_MAKE_MULTI_SPECIES_DIFFUSION_ENGINE_H
#define DIFFUSION3D_MAKE_MULTI_SPECIES_DIFFUSION_ENGINE_H

#include <memory>
#include "multi_species_diffusion_engine.h"
#include "explicit_multi_species_heat_stepper.h"

/**
 * @brief Factory function to create a MultiSpeciesDiffusionEngine.
 * @param species List of species IDs.
 * @param settings Diffusion settings.
 * @return Unique pointer to a MultiSpeciesDiffusionEngine.
 */
inline std::unique_ptr<MultiSpeciesDiffusionEngine> MakeMultiSpeciesDiffusionEngine(
    const std::vector<SpeciesId> &species,
    const MultiSpeciesDiffusionSettings &settings)
{
    auto engine = std::make_unique<ExplicitMultiSpeciesHeatStepper>();
    engine->configure_species_interval(MultiSpeciesFieldGrid(species, 1, 1, 1), settings); // Dummy grid for validation
    return engine;
}

#endif // DIFFUSION3D_MAKE_MULTI_SPECIES_DIFFUSION_ENGINE_H
