/**
 * @file scalar_field_grid.h
 * @brief 3D scalar field grid for diffusion in Diffusion3D (adapted from Diffusion3DContext).
 *
 * Encapsulates a single 3D double-precision grid with row-major C layout.
 * Replaces the (u, u_next) dual-buffer model from Diffusion3DContext.
 */

#ifndef DIFFUSION3D_SCALAR_FIELD_GRID_H
#define DIFFUSION3D_SCALAR_FIELD_GRID_H

#include <vector>
#include <cassert>

/**
 * @class ScalarFieldGrid
 * @brief Single species' 3D scalar field on uniform Cartesian grid.
 *
 * Row-major C layout: data[i] = u(x,y,z) where i = x + nx*(y + ny*z).
 * Matches the 6-point stencil kernels from diffusion3d.
 */
class ScalarFieldGrid
{
public:
    int nx_, ny_, nz_;         ///< Grid dimensions (same units as h)
    std::vector<double> data_; ///< Flattened grid: u[x + nx*(y + ny*z)]

    /**
     * @brief Construct a grid with given dimensions, initialized to zero.
     * @param nx, ny, nz Grid dimensions
     */
    ScalarFieldGrid(int nx, int ny, int nz)
        : nx_(nx), ny_(ny), nz_(nz), data_(static_cast<std::size_t>(nx) * ny * nz, 0.0)
    {
        assert(nx > 0 && ny > 0 && nz > 0);
    }

    /**
     * @brief Linear index for (x, y, z) in row-major order.
     */
    inline int idx(int x, int y, int z) const
    {
        assert(x >= 0 && x < nx_ && y >= 0 && y < ny_ && z >= 0 && z < nz_);
        return x + nx_ * (y + ny_ * z);
    }

    /**
     * @brief Total number of grid points.
     */
    inline int size() const { return nx_ * ny_ * nz_; }

    /**
     * @brief Mutable access to element at (x, y, z).
     */
    double &at(int x, int y, int z)
    {
        return data_[idx(x, y, z)];
    }

    /**
     * @brief Const access to element at (x, y, z).
     */
    const double &at(int x, int y, int z) const
    {
        return data_[idx(x, y, z)];
    }
};

#endif // DIFFUSION3D_SCALAR_FIELD_GRID_H
