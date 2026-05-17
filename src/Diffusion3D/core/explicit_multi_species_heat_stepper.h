#ifndef DIFFUSION3D_EXPLICIT_MULTI_SPECIES_HEAT_STEPPER_H
#define DIFFUSION3D_EXPLICIT_MULTI_SPECIES_HEAT_STEPPER_H

// Explicit multi-species heat equation stepper
#include "multi_species_diffusion_engine.h"

class ExplicitMultiSpeciesHeatStepper : public MultiSpeciesDiffusionEngine
{
public:
    /**
     * @brief Configure the stepper for a given species interval.
     */
    void configure_species_interval(const MultiSpeciesFieldGrid &grid, const MultiSpeciesDiffusionSettings &settings) override;
    /**
     * @brief Advance all species by one interval.
     */
    void advance_species_interval(MultiSpeciesFieldGrid &grid) override;
    /**
     * @brief Get the supported species IDs.
     */
    std::vector<SpeciesId> supported_species() const override;
    /**
     * @brief Get the resolved algorithm type.
     */
    DiffusionAlgorithm resolved_algorithm() const override { return DiffusionAlgorithm::ExplicitHeatEquation; }

private:
    MultiSpeciesDiffusionSettings settings_;
    std::vector<SpeciesId> species_ids_;
};

#endif // DIFFUSION3D_EXPLICIT_MULTI_SPECIES_HEAT_STEPPER_H
