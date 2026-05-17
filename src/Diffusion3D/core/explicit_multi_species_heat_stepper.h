#ifndef DIFFUSION3D_EXPLICIT_MULTI_SPECIES_HEAT_STEPPER_H
#define DIFFUSION3D_EXPLICIT_MULTI_SPECIES_HEAT_STEPPER_H

/**
 * @file explicit_multi_species_heat_stepper.h
 * @brief CPU explicit Euler implementation of MultiSpeciesDiffusionEngine (adapted from DiffusionSolver).
 *
 * Implements multi-species diffusion PDE stepping using explicit Euler with 6-point stencil.
 * Handles per-species CFL substep planning for mixed diffusivities.
 */

#include <vector>
#include <map>
#include "multi_species_diffusion_engine.h"

/**
 * @class ExplicitMultiSpeciesHeatStepper
 * @brief CPU explicit Euler stepper for multi-species heat equation (Phase I reference implementation).
 *
 * **Algorithm**: Standard explicit Euler with 6-point stencil:
 *   u_new[i] = u[i] + dt_sub * D/h^2 * (sum_neighbors - 6*u[i])
 *
 * **CFL Stability**: dt_max = (h^2 / 6D) * safety
 * For mixed diffusivities, uses dt_min across all species (conservative).
 *
 * **Per-Species Substeps**: Plans substep count N_i = ceil(Δt / dt_i,max) for species i.
 * Executes N_i substeps at dt = Δt/N_i for each species independently, ensuring all stay stable.
 *
 * **Boundary Conditions**: Dirichlet (u=0 at domain walls).
 * Interior-point-only kernels; boundary remains zero throughout stepping.
 */
class ExplicitMultiSpeciesHeatStepper : public MultiSpeciesDiffusionEngine
{
public:
    ExplicitMultiSpeciesHeatStepper() = default;
    virtual ~ExplicitMultiSpeciesHeatStepper() = default;

    /**
     * @brief Configure the stepper for a diffusion step (same as DiffusionSolver::configure_tick).
     *
     * Precomputes CFL substep plan (n_sub, dt_sub) and validates grid/settings consistency.
     * Must be called once per macro timestep before advance_species_interval().
     *
     * @param grid MultiSpeciesFieldGrid with configured species and dimensions
     * @param settings Diffusion parameters (species_diffusivities, safety, etc.)
     *
     * @throws std::runtime_error if settings invalid or species mismatch
     */
    void configure_species_interval(
        const MultiSpeciesFieldGrid &grid,
        const MultiSpeciesDiffusionSettings &settings) override;

    /**
     * @brief Advance all species by the configured global timestep.
     *
     * For each species: execute N_i substeps at dt_sub = dt_global / N_i.
     * All interior points updated via 6-point stencil; boundary remains zero.
     *
     * @param grid MultiSpeciesFieldGrid (modified in-place)
     *
     * @throws std::runtime_error if not configured before call
     */
    void advance_species_interval(MultiSpeciesFieldGrid &grid) override;

    /**
     * @brief Get configured species IDs.
     * @return Vector from last configure_species_interval() call
     */
    std::vector<SpeciesId> supported_species() const override;

    /**
     * @brief Get algorithm type.
     * @return Always DiffusionAlgorithm::ExplicitHeatEquation
     */
    DiffusionAlgorithm resolved_algorithm() const override
    {
        return DiffusionAlgorithm::ExplicitHeatEquation;
    }

    // Accessors for testing/debugging
    double tick_dt() const { return tick_dt_; }
    int n_sub() const { return n_sub_; }
    double dt_sub() const { return dt_sub_; }

private:
    /**
     * @brief Apply 6-point stencil kernel to one species (adapted from cpu_diffusion3d.cpp).
     *
     * @param grid ScalarFieldGrid to update (modified in-place)
     * @param D Diffusivity coefficient
     * @param h Grid spacing
     * @param dt_sub Micro-timestep for this stencil application
     */
    void apply_stencil_kernel(
        ScalarFieldGrid &grid,
        double D,
        double h,
        double dt_sub);

    // Configuration state
    bool configured_ = false;
    double tick_dt_ = 0.0;
    int n_sub_ = 1;
    double dt_sub_ = 0.0;

    // Grid parameters (validated on configure)
    int cfg_nx_ = 0, cfg_ny_ = 0, cfg_nz_ = 0;
    double cfg_h_ = 1.0;

    // Configured settings and species list
    MultiSpeciesDiffusionSettings settings_;
    std::vector<SpeciesId> species_ids_;

    // Pre-computed dt_max for each species (for validation only)
    std::map<SpeciesId, double> dt_max_per_species_;
};

#endif // DIFFUSION3D_EXPLICIT_MULTI_SPECIES_HEAT_STEPPER_H
