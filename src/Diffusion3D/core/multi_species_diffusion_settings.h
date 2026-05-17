#ifndef DIFFUSION3D_MULTI_SPECIES_DIFFUSION_SETTINGS_H
#define DIFFUSION3D_MULTI_SPECIES_DIFFUSION_SETTINGS_H

#include <cassert>
#include <map>
#include <vector>
#include <stdexcept>

using SpeciesId = int;

/**
 * @file multi_species_diffusion_settings.h
 * @brief Multi-species diffusion configuration (adapted from DiffusionParams).
 *
 * Stores per-species diffusivity, safety factor, and backend selection.
 * Each species has independent diffusivity D with unified safety scaling.
 */

/**
 * @enum DiffusionBackend
 * @brief Hardware/algorithm backend selection (same as diffusion3d_solver.h).
 */
enum class DiffusionBackend
{
    CpuStencil,
    GpuStencil,
    GpuFftPrecomputed,
    Auto,
};

/**
 * @brief Multi-species diffusion parameters (adapted from DiffusionParams).
 *
 * Stores per-species diffusivity D, unified safety factor, backend choice.
 * All species share the same physical domain (nx, ny, nz) and spacing h.
 */
struct MultiSpeciesDiffusionSettings
{
    /**
     * @brief Per-species diffusivity D in du_i/dt = D_i * Laplacian(u_i).
     *
     * Maps SpeciesId -> diffusivity value. All D > 0 required.
     * Mixed diffusivities constrain global timestep: dt_max = min_i((h^2 / (6*D_i)) * safety).
     */
    std::map<SpeciesId, double> species_diffusivities;

    /**
     * @brief Optional species names (for logging, debugging).
     */
    std::vector<std::string> species_names;

    /**
     * @brief CFL safety factor: 0 < safety <= 1.0
     *
     * Scales the explicit-Euler stability bound.
     * Lower = more conservative (smaller steps, more substeps).
     */
    double safety = 1.0;

    /**
     * @brief Hardware backend selection (same enum as diffusion3d).
     */
    DiffusionBackend backend = DiffusionBackend::Auto;

    /**
     * @brief FFT padding extents (same as DiffusionParams, for GPU FFT backend).
     */
    int fft_real_extent_x = 0;
    int fft_real_extent_y = 0;
    int fft_real_extent_z = 0;

    /**
     * @brief Validate all parameters for consistency.
     *
     * @throws std::runtime_error if:
     *   - species_diffusivities is empty
     *   - any diffusivity <= 0
     *   - safety not in (0, 1]
     */
    void validate() const
    {
        if (species_diffusivities.empty())
            throw std::runtime_error("species_diffusivities is empty");

        for (const auto &[species_id, D] : species_diffusivities)
        {
            assert(species_id >= 0);
            if (D <= 0.0)
                throw std::runtime_error("species_diffusivities[" + std::to_string(species_id) +
                                         "] = " + std::to_string(D) + " but must be > 0");
        }

        if (safety <= 0.0 || safety > 1.0)
            throw std::runtime_error("safety = " + std::to_string(safety) + " but must be in (0, 1]");
    }
};

#endif // DIFFUSION3D_MULTI_SPECIES_DIFFUSION_SETTINGS_H
