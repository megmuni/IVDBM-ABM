/**
 * @file explicit_multi_species_heat_stepper.cpp
 * @brief CPU explicit Euler multi-species diffusion.
 */

#include "explicit_multi_species_heat_stepper.h"
#include "diffusion3d_step_euler_cpu.h"
#include <stdexcept>

void ExplicitMultiSpeciesHeatStepper::configure_species_interval(
    const MultiSpeciesFieldGrid &grid,
    const MultiSpeciesDiffusionSettings &settings,
    double tick_dt)
{
    settings.validate();

    species_ids_ = grid.species();
    if (species_ids_.empty())
        throw std::runtime_error("Grid must contain at least one species");

    for (auto id : species_ids_)
    {
        if (settings.species_diffusivities.find(id) == settings.species_diffusivities.end())
            throw std::runtime_error("Settings missing species " + std::to_string(id));
    }

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

    configured_ = true;
}

void ExplicitMultiSpeciesHeatStepper::advance_species_interval(MultiSpeciesFieldGrid &grid)
{
    if (!configured_)
        throw std::runtime_error("Stepper not configured; call configure_species_interval() first");

    if (grid.nx_ != cfg_nx_ || grid.ny_ != cfg_ny_ || grid.nz_ != cfg_nz_ || grid.h_ != cfg_h_)
        throw std::runtime_error("Grid dimensions changed since configure_species_interval()");

    for (auto species_id : species_ids_)
    {
        ScalarFieldGrid &species_grid = grid.grid(species_id);
        const double D = settings_.species_diffusivities.at(species_id);
        const SubstepPlan &plan = plan_per_species_.at(species_id);

        for (int i = 0; i < plan.n_sub; ++i)
            diffusion3d_step_euler_scalar(species_grid, cfg_h_, D, plan.dt_sub);
    }
}

std::vector<SpeciesId> ExplicitMultiSpeciesHeatStepper::supported_species() const
{
    return species_ids_;
}

int ExplicitMultiSpeciesHeatStepper::n_sub(SpeciesId id) const
{
    return plan_per_species_.at(id).n_sub;
}

double ExplicitMultiSpeciesHeatStepper::dt_sub(SpeciesId id) const
{
    return plan_per_species_.at(id).dt_sub;
}
