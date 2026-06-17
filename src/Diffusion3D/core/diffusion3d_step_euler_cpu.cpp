/**
 * @file diffusion3d_step_euler_cpu.cpp
 * @brief CPU reference: explicit Euler plus 6-point face Laplacian (open box).
 */

#include "diffusion3d_step_euler_cpu.h"

#include <cassert>

namespace {

inline bool is_inside_domain(int x, int y, int z, int nx, int ny, int nz)
{
    return x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz;
}

void apply_open_box_euler_step(const std::vector<double> &u,
                               std::vector<double> &u_next,
                               int nx, int ny, int nz,
                               double h, double D, double dt_sub)
{
    assert(h > 0.0);
    assert(D >= 0.0);
    assert(dt_sub >= 0.0);
    if (dt_sub == 0.0)
        return;

    const double lap_scale = D / (h * h);

    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const int center_idx = x + nx * (y + ny * z);
                const double center = u[static_cast<size_t>(center_idx)];
                double lap = 0.0;

                const int nbr[6][3] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                                       {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};

                for (int k = 0; k < 6; ++k)
                {
                    const int nx_ = x + nbr[k][0];
                    const int ny_ = y + nbr[k][1];
                    const int nz_ = z + nbr[k][2];

                    if (is_inside_domain(nx_, ny_, nz_, nx, ny, nz))
                    {
                        const int nidx = nx_ + nx * (ny_ + ny * nz_);
                        lap += u[static_cast<size_t>(nidx)] - center;
                    }
                }

                u_next[static_cast<size_t>(center_idx)] = center + lap_scale * dt_sub * lap;
            }
        }
    }
}

} // namespace

void diffusion3d_step_euler_scalar(ScalarFieldGrid &grid, double h, double D, double dt_sub)
{
    if (dt_sub == 0.0)
        return;
    grid.ensure_scratch();
    apply_open_box_euler_step(grid.data_, grid.scratch_, grid.nx_, grid.ny_, grid.nz_, h, D, dt_sub);
    grid.data_.swap(grid.scratch_);
}
