#ifndef DIFFUSION3D_MULTI_SPECIES_DIFFUSION_ENGINE_H
#define DIFFUSION3D_MULTI_SPECIES_DIFFUSION_ENGINE_H

#include <vector>
#include "multi_species_field_grid.h"
#include "multi_species_diffusion_settings.h"

using SpeciesId = int;

/**
 * @file multi_species_diffusion_engine.h
 * @brief Abstract interface for multi-species diffusion PDE stepping (adapted from DiffusionSolver).
 *
 * Defines MultiSpeciesDiffusionEngine (Strategy pattern) and DiffusionAlgorithm enum.
 * Enables CPU explicit, GPU stencil, GPU FFT backends with identical calling convention.
 */

/**
 * @enum DiffusionAlgorithm
 * @brief PDE stepping algorithm (same enum values as diffusion3d_solver.h DiffusionBackend).
 */
enum class DiffusionAlgorithm
{
    ExplicitHeatEquation = 0,  ///< CPU explicit Euler with 6-point stencil
    GpuStencil = 1,             ///< GPU 6-point stencil (Phase II)
    GpuFftPrecomputed = 2,      ///< GPU spectral via FFT (Phase III)
};

/**
 * @class MultiSpeciesDiffusionEngine
 * @brief Abstract base for multi-species diffusion PDE solvers (Strategy pattern).
 *
 * Encapsulates algorithm choice (CPU explicit, GPU stencil, GPU FFT) and CFL planning.
 * Each species has independent diffusivity D but shares global timestep dt_global.
 */
class MultiSpeciesDiffusionEngine
{
public:
    virtual ~MultiSpeciesDiffusionEngine() = default;

    /**
     * @brief Configure the engine for a diffusion step.
     *
     * Precomputes CFL timesteps, allocates GPU buffers (if needed), validates settings.
     * Must be called once per tick before advance_species_interval().
     *
     * @param grid MultiSpeciesFieldGrid with configured species
     * @param settings Multi-species diffusivities and safety factor
     * @throws std::runtime_error on configuration failure
     */
    virtual void configure_species_interval(
        const MultiSpeciesFieldGrid &grid,
        const MultiSpeciesDiffusionSettings &settings) = 0;

    /**
     * @brief Advance all species by their configured macro timestep.
     *
     * Executes per-species sub-stepping (CFL-stable dt_sub) for all species.
     * Interior points updated via stencil kernel; boundary remains zero (Dirichlet).
     *
     * @param grid MultiSpeciesFieldGrid (modified in-place)
     * @throws std::runtime_error on execution failure (GPU errors, etc.)
     */
    virtual void advance_species_interval(MultiSpeciesFieldGrid &grid) = 0;

    /**
     * @brief Get the species IDs configured in the last configure_species_interval().
     * @return Vector of species_ids from grid
     */
    virtual std::vector<SpeciesId> supported_species() const = 0;

    /**
     * @brief Get the resolved algorithm type.
     * @return DiffusionAlgorithm enum (never Auto)
     */
    virtual DiffusionAlgorithm resolved_algorithm() const = 0;
};

#endif // DIFFUSION3D_MULTI_SPECIES_DIFFUSION_ENGINE_H
