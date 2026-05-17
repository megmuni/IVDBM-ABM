#include "diffusion3d_common.h"

#include <cassert>

/**
 * @file cpu_diffusion3d.cpp
 * @brief CPU reference: explicit Euler plus 6-point face Laplacian (open box: missing neighbors omitted).
 *
 * PDE: du/dt = D * Laplacian(u) on a uniform grid with spacing h.
 * Laplacian: sum over six face neighbors of (u_neighbor - u_center) / h^2 (centered second difference).
 * Update: u_new = u + dt_sub * D/h^2 * sum_face (u_neighbor - u_center).
 */

/**
 * @brief True if (x,y,z) lies inside [0,nx) x [0,ny) x [0,nz).
 */
static inline bool is_inside_domain(int x, int y, int z, int nx, int ny, int nz)
{
    return x >= 0 && x < nx &&
           y >= 0 && y < ny &&
           z >= 0 && z < nz;
}

void diffusion3d_step_euler_cpu(Diffusion3DContext &ctx, double dt_sub)
{
    assert(ctx.h > 0.0);
    assert(ctx.D >= 0.0);
    assert(dt_sub >= 0.0);

    const double inv_h2 = 1.0 / (ctx.h * ctx.h);
    const double lap_scale = ctx.D * inv_h2; // D / h^2

    for (int z = 0; z < ctx.nz; ++z)
    {
        for (int y = 0; y < ctx.ny; ++y)
        {
            for (int x = 0; x < ctx.nx; ++x)
            {
                const double center = ctx.u[ctx.idx(x, y, z)];
                double lap = 0.0;

                const int nbr[6][3] = {
                    { 1, 0, 0}, {-1, 0, 0},
                    { 0, 1, 0}, { 0,-1, 0},
                    { 0, 0, 1}, { 0, 0,-1}
                };

                for (int k = 0; k < 6; ++k)
                {
                    const int nx_ = x + nbr[k][0];
                    const int ny_ = y + nbr[k][1];
                    const int nz_ = z + nbr[k][2];

                    if (is_inside_domain(nx_, ny_, nz_, ctx.nx, ctx.ny, ctx.nz))
                    {
                        const int nidx = ctx.idx(nx_, ny_, nz_);
                        lap += ctx.u[nidx] - center;
                    }
                }

                ctx.u_next[ctx.idx(x, y, z)] = center + lap_scale * dt_sub * lap;
            }
        }
    }

    ctx.u.swap(ctx.u_next);
}
