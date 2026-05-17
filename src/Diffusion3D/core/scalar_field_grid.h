// ScalarFieldGrid: 3D double array for diffusion
#include <vector>

/**
 * @file scalar_field_grid.h
 * @brief 3D scalar field grid for diffusion in Diffusion3D.
 *
 * Defines ScalarFieldGrid for 3D double-precision arrays.
 */

#ifndef DIFFUSION3D_SCALAR_FIELD_GRID_H
#define DIFFUSION3D_SCALAR_FIELD_GRID_H

#include <vector>

/**
 * @class ScalarFieldGrid
 * @brief 3D double array for scalar field diffusion.
 */
class ScalarFieldGrid
{
public:
    int nx_, ny_, nz_;
    std::vector<double> data_;

    ScalarFieldGrid(int nx, int ny, int nz)
        : nx_(nx), ny_(ny), nz_(nz), data_(nx * ny * nz, 0.0) {}

    double &at(int x, int y, int z)
    {
        return data_[x + nx_ * (y + ny_ * z)];
    }

    const double &at(int x, int y, int z) const
    {
        return data_[x + nx_ * (y + ny_ * z)];
    }
};

#endif // DIFFUSION3D_SCALAR_FIELD_GRID_H
