/**
 * @file explicit_multi_species_heat_stepper.cpp
 * @brief CPU explicit Euler implementation for multi-species diffusion (adapted from diffusion3d).
 *
 * Implements ExplicitMultiSpeciesHeatStepper: CFL planning, per-species substep execution,
 * and 6-point stencil kernel (adapted from cpu_diffusion3d.cpp).
 */

#include "explicit_multi_species_heat_stepper.h"
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>

namespace
{

/**
 * @brief Compute CFL stability bound: dt_max = safety * h^2 / (6*D).
 *
 * Adapted from diffusion3d_timestep.h::compute_stability_constraint.
 */
double compute_stability_constraint(double h, double D, double safety)
{
    assert(h > 0.0 && D > 0.0 && safety > 0.0 && safety <= 1.0);
    return safety * (h * h) / (6.0 * D);
}

/**
 * @brief Plan substeps: n_sub * dt_sub = tick_dt, all dt_sub <= dt_max.
 *
 * Adapted from diffusion3d_timestep.h::plan_substeps.
 */
struct SubstepPlan
{
    int n_sub;
    double dt_sub;
};

SubstepPlan plan_substeps(double tick_dt, double dt_max)
{
    assert(tick_dt >= 0.0);
    if (tick_dt == 0.0)
        return {1, 0.0};
    if (!std::isfinite(dt_max))
        return {1, tick_dt};
    assert(dt_max > 0.0);

    const double ratio = tick_dt / dt_max;
    const int n_sub = std::max(1, static_cast<int>(std::ceil(ratio)));
    const double dt_sub = tick_dt / static_cast<double>(n_sub);
    return {n_sub, dt_sub};
}

/**
 * @brief Check if (x,y,z) is inside [0,nx) x [0,ny) x [0,nz).
 */
inline bool is_inside_domain(int x, int y, int z, int nx, int ny, int nz)
{
    return x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz;
}

} // namespace

void ExplicitMultiSpeciesHeatStepper::configure_species_interval(
    const MultiSpeciesFieldGrid &grid,
    const MultiSpeciesDiffusionSettings &settings)
{
    settings.validate();

    species_ids_ = grid.species();
    if (species_ids_.empty())
        throw std::runtime_error("Grid must contain at least one species");

    // Ensure all species in grid have settings
    for (auto id : species_ids_)
    {
        if (settings.species_diffusivities.find(id) == settings.species_diffusivities.end())
            throw std::runtime_error("Settings missing species " + std::to_string(id));
    }

    // Store configuration
    settings_ = settings;
    cfg_nx_ = grid.nx_;
    cfg_ny_ = grid.ny_;
    cfg_nz_ = grid.nz_;
    cfg_h_ = grid.h_;

    // Compute CFL bounds for each species
    double dt_max_global = std::numeric_limits<double>::infinity();
    dt_max_per_species_.clear();

    for (auto id : species_ids_)
    {
        double D = settings.species_diffusivities.at(id);
        double dt_max = compute_stability_constraint(cfg_h_, D, settings.safety);
        dt_max_per_species_[id] = dt_max;
        dt_max_global = std::min(dt_max_global, dt_max);
    }

    // For Phase I, use uniform global timestep (conservative)
    // dt_global = min_i(dt_max_i)
    // This ensures all species stay stable with single substep plan

    configured_ = true;
}

void ExplicitMultiSpeciesHeatStepper::advance_species_interval(MultiSpeciesFieldGrid &grid)
{
    if (!configured_)
        throw std::runtime_error("Stepper not configured; call configure_species_interval() first");

    // Validate grid consistency
    if (grid.nx_ != cfg_nx_ || grid.ny_ != cfg_ny_ || grid.nz_ != cfg_nz_ || grid.h_ != cfg_h_)
        throw std::runtime_error("Grid dimensions changed since configure_species_interval()");

    // Compute global timestep: dt_global = min_i(dt_max_i)
    // This is conservative but simplifies Phase I
    double dt_global = std::numeric_limits<double>::infinity();
    for (const auto &[id, dt_max] : dt_max_per_species_)
        dt_global = std::min(dt_global, dt_max);

    if (!std::isfinite(dt_global) || dt_global <= 0.0)
        return;  // No diffusion or zero timestep

    // Advance each species with per-species substep planning
    for (auto species_id : species_ids_)
    {
        ScalarFieldGrid &species_grid = grid.grid(species_id);
        double D = settings_.species_diffusivities.at(species_id);
        double dt_max = dt_max_per_species_.at(species_id);

        // Plan substeps for this species to stay CFL stable
        SubstepPlan plan = plan_substeps(dt_global, dt_max);
        n_sub_ = plan.n_sub;
        dt_sub_ = plan.dt_sub;
        tick_dt_ = dt_global;

        // Execute substeps
        for (int i = 0; i < n_sub_; ++i)
            apply_stencil_kernel(species_grid, D, cfg_h_, dt_sub_);
    }
}

std::vector<SpeciesId> ExplicitMultiSpeciesHeatStepper::supported_species() const
{
    return species_ids_;
}

void ExplicitMultiSpeciesHeatStepper::apply_stencil_kernel(
    ScalarFieldGrid &grid,
    double D,
    double h,
    double dt_sub)
{
    // Validate
    assert(h > 0.0 && D >= 0.0 && dt_sub >= 0.0);
    if (dt_sub == 0.0)
        return;

    const int nx = grid.nx_;
    const int ny = grid.ny_;
    const int nz = grid.nz_;
    const double inv_h2 = 1.0 / (h * h);
    const double coeff = D * inv_h2 * dt_sub;  // D * dt_sub / h^2

    // Temporary buffer for u_new
    std::vector<double> u_next = grid.data_;

    // Apply 6-point stencil to interior points (adapted from cpu_diffusion3d.cpp)
    // Boundary points (walls) remain zero (Dirichlet)
    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const int idx = grid.idx(x, y, z);
                const double center = grid.data_[idx];
                double laplacian = 0.0;

                // Six neighbors: x±1, y±1, z±1
                const int nbr[6][3] = {
                    {1, 0, 0},
                    {-1, 0, 0},
                    {0, 1, 0},
                    {0, -1, 0},
                    {0, 0, 1},
                    {0, 0, -1}};

                for (int k = 0; k < 6; ++k)
                {
                    const int nx_ = x + nbr[k][0];
                    const int ny_ = y + nbr[k][1];
                    const int nz_ = z + nbr[k][2];

                    if (is_inside_domain(nx_, ny_, nz_, nx, ny, nz))
                    {
                        const int nidx = grid.idx(nx_, ny_, nz_);
                        laplacian += grid.data_[nidx] - center;
                    }
                }

                // Explicit Euler: u_new = u + coeff * laplacian
                u_next[idx] = center + coeff * laplacian;
            }
        }
    }

    // Swap buffers
    grid.data_ = u_next;
}
