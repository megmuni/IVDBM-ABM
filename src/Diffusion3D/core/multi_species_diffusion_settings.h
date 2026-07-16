#ifndef DIFFUSION3D_MULTI_SPECIES_DIFFUSION_SETTINGS_H
#define DIFFUSION3D_MULTI_SPECIES_DIFFUSION_SETTINGS_H

#include "diffusion_algorithm.h"
#include <cassert>
#include <map>
#include <stdexcept>
#include <vector>

using SpeciesId = int;

/**
 * @file multi_species_diffusion_settings.h
 * @brief Multi-species diffusion configuration.
 */

/**
 * @brief Multi-species diffusion parameters.
 *
 * Stores per-species diffusivity D, unified safety factor, backend choice.
 * All species share the same physical domain (nx, ny, nz) and spacing h.
 */
struct MultiSpeciesDiffusionSettings {
  /**
   * @brief Per-species diffusivity D in du_i/dt = D_i * Laplacian(u_i).
   *
   * Maps SpeciesId -> diffusivity value. All D > 0 required.
   * Mixed diffusivities constrain global timestep: dt_max = min_i((h^2 /
   * (6*D_i)) * safety).
   */
  std::map<SpeciesId, double> species_diffusivities;

  /**
   * @brief Optional species names (for logging, debugging).
   */
  std::vector<std::string> species_names;

  /**
   * @brief CFL safety factor: 0 < safety < 1.0
   *
   * Scales the explicit-Euler stability bound.
   * Lower = more conservative (smaller steps, more substeps).
   *
   * Must be strictly < 1. At safety == 1 the substep lands exactly on the
   * stability edge (s = D*dt_sub/h^2 = 1/6). Requiring safety < 1 keeps
   * the centre weight >= 1 - safety > 0 for every D, h and tick.
   */
  double safety = 0.5;

  /**
   * @brief Algorithm / hardware backend selection.
   */
  DiffusionAlgorithm algorithm = DiffusionAlgorithm::ExplicitHeatEquation;

  /**
   * @brief Validate all parameters for consistency.
   *
   * @throws std::runtime_error if:
   *   - species_diffusivities is empty
   *   - any diffusivity <= 0
   *   - safety not in (0, 1)
   */
  void validate() const {
    if (species_diffusivities.empty())
      throw std::runtime_error("species_diffusivities is empty");

    for (const auto &[species_id, D] : species_diffusivities) {
      assert(species_id >= 0);
      if (D <= 0.0)
        throw std::runtime_error(
            "species_diffusivities[" + std::to_string(species_id) +
            "] = " + std::to_string(D) + " but must be > 0");
    }

    if (safety <= 0.0 || safety >= 1.0)
      throw std::runtime_error("safety = " + std::to_string(safety) +
                               " but must be in (0, 1)");
  }
};

#endif // DIFFUSION3D_MULTI_SPECIES_DIFFUSION_SETTINGS_H
