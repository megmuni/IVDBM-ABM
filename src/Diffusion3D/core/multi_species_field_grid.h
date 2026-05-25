/**
 * @file multi_species_field_grid.h
 * @brief Multi-species grid container on a shared spatial domain.
 */

#ifndef DIFFUSION3D_MULTI_SPECIES_FIELD_GRID_H
#define DIFFUSION3D_MULTI_SPECIES_FIELD_GRID_H

#include <map>
#include <memory>
#include <vector>
#include <cassert>
#include "scalar_field_grid.h"

using SpeciesId = int;

/**
 * @class MultiSpeciesFieldGrid
 * @brief Container for all species' 3D scalar fields on shared domain (nx, ny, nz).
 *
 * Each species owns a ScalarFieldGrid; all share dimensions and spacing h.
 */
class MultiSpeciesFieldGrid
{
public:
    int nx_, ny_, nz_;                                                    ///< Shared grid dimensions
    double h_;                                                            ///< Shared grid spacing
    std::map<SpeciesId, std::shared_ptr<ScalarFieldGrid>> species_grids_; ///< Per-species grids
    std::vector<SpeciesId> species_ids_;                                  ///< Species list (order)

    /**
     * @brief Construct a multi-species grid with shared dimensions.
     * @param ids List of species IDs (must be non-empty and unique)
     * @param nx, ny, nz Grid dimensions
     * @param h Grid spacing
     */
    MultiSpeciesFieldGrid(const std::vector<SpeciesId> &ids, int nx, int ny, int nz, double h = 1.0)
        : nx_(nx), ny_(ny), nz_(nz), h_(h), species_ids_(ids)
    {
        assert(!ids.empty());
        assert(nx > 0 && ny > 0 && nz > 0 && h > 0);
        for (auto id : ids)
        {
            species_grids_[id] = std::make_shared<ScalarFieldGrid>(nx, ny, nz);
        }
    }

    /**
     * @brief Mutable access to a species' grid.
     * @throws std::out_of_range if species not in grid
     */
    ScalarFieldGrid &grid(SpeciesId id)
    {
        return *species_grids_.at(id);
    }

    /**
     * @brief Const access to a species' grid.
     * @throws std::out_of_range if species not in grid
     */
    const ScalarFieldGrid &grid(SpeciesId id) const
    {
        return *species_grids_.at(id);
    }

    /**
     * @brief Check if species exists in grid.
     */
    bool has_species(SpeciesId id) const
    {
        return species_grids_.count(id) > 0;
    }

    /**
     * @brief Get list of all species.
     */
    const std::vector<SpeciesId> &species() const
    {
        return species_ids_;
    }
};

#endif // DIFFUSION3D_MULTI_SPECIES_FIELD_GRID_H
