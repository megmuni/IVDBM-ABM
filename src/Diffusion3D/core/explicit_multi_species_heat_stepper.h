#ifndef DIFFUSION3D_EXPLICIT_MULTI_SPECIES_HEAT_STEPPER_H
#define DIFFUSION3D_EXPLICIT_MULTI_SPECIES_HEAT_STEPPER_H

/**
 * @file explicit_multi_species_heat_stepper.h
 * @brief CPU explicit Euler stepper for multi-species heat equation.
 */

#include <vector>
#include <map>
#include "multi_species_diffusion_engine.h"
#include "diffusion3d_timestep.h"

/**
 * @class ExplicitMultiSpeciesHeatStepper
 * @brief CPU explicit Euler stepper for multi-species heat equation (Phase I reference).
 *
 * Matches configure_species_interval + advance_species_interval with per-species substep plans.
 */
class ExplicitMultiSpeciesHeatStepper : public MultiSpeciesDiffusionEngine
{
public:
    ExplicitMultiSpeciesHeatStepper() = default;
    virtual ~ExplicitMultiSpeciesHeatStepper() = default;

    void configure_species_interval(
        const MultiSpeciesFieldGrid &grid,
        const MultiSpeciesDiffusionSettings &settings,
        double tick_dt) override;

    void advance_species_interval(MultiSpeciesFieldGrid &grid) override;

    std::vector<SpeciesId> supported_species() const override;

    DiffusionAlgorithm resolved_algorithm() const override
    {
        return DiffusionAlgorithm::ExplicitHeatEquation;
    }

    double tick_dt() const { return tick_dt_; }
    int n_sub(SpeciesId id) const;
    double dt_sub(SpeciesId id) const;

private:
    bool configured_ = false;
    double tick_dt_ = 0.0;

    int cfg_nx_ = 0, cfg_ny_ = 0, cfg_nz_ = 0;
    double cfg_h_ = 1.0;

    MultiSpeciesDiffusionSettings settings_;
    std::vector<SpeciesId> species_ids_;
    std::map<SpeciesId, SubstepPlan> plan_per_species_;
};

#endif // DIFFUSION3D_EXPLICIT_MULTI_SPECIES_HEAT_STEPPER_H
