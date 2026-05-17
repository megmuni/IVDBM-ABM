/**
 * @file multi_species_field_grid.h
 * @brief Multi-species grid container with memory mapping for Diffusion3D.
 *
 * Defines MultiSpeciesFieldGrid for managing per-species scalar fields.
 */

#ifndef DIFFUSION3D_MULTI_SPECIES_FIELD_GRID_H
#define DIFFUSION3D_MULTI_SPECIES_FIELD_GRID_H

#include <map>
#include <memory>
#include <vector>
#include "scalar_field_grid.h"

using SpeciesId = int; // Placeholder, replace with canonical SpeciesId later

/**
 * @class MultiSpeciesFieldGrid
 * @brief Container for multiple species' scalar field grids.
 */
class MultiSpeciesFieldGrid
{
public:
    std::map<SpeciesId, std::shared_ptr<ScalarFieldGrid>> species_grids_;
    std::vector<SpeciesId> species_ids_;

    /**
     * @brief Construct a grid for each species.
     */
    MultiSpeciesFieldGrid(const std::vector<SpeciesId> &ids, int nx, int ny, int nz) : species_ids_(ids)
    {
        for (auto id : ids)
        {
            species_grids_[id] = std::make_shared<ScalarFieldGrid>(nx, ny, nz);
        }
    }

    /**
     * @brief Access a species' grid.
     */
    ScalarFieldGrid &grid(SpeciesId id)
    {
        return *species_grids_.at(id);
    }

    /**
     * @brief Const access to a species' grid.
     */
    const ScalarFieldGrid &grid(SpeciesId id) const
    {
        return *species_grids_.at(id);
    }
};

#endif // DIFFUSION3D_MULTI_SPECIES_FIELD_GRID_H
