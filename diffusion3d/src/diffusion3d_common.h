#ifndef DIFFUSION3D_COMMON_H
#define DIFFUSION3D_COMMON_H

#include <cassert>
#include <vector>

/**
 * @file diffusion3d_common.h
 * @brief Shared context and API for 3D diffusion (CPU & GPU implementations).
 *
 * This file defines the data structures used by both CPU and GPU solvers
 * for simulating the 3D heat equation using a 6-point finite-difference stencil.
 *
 * Prefer `Diffusion3DContext::make(nx, ny, nz, h, D)` for new code; use
 * `make_with_reference_dt` only when you want a demo/reference timestep stored on the context.
 * `DiffusionSolver` does not read `stored_dt`.
 *
 * References:
 * 1. https://student.cs.uwaterloo.ca/~cs475/CS475-Lecture04.pdf
 * 2. https://enccs.github.io/OpenACC-CUDA-beginners/2.02_cuda-heat-equation/
 */

/**
 * @brief Simulation state for a 3D diffusion (heat equation) solver.
 *
 * This structure stores:
 * - The spatial grid dimensions
 * - GPU execution configuration (block size)
 * - The scalar fields: current state `u` and Euler scratch buffer `u_next`
 * - Physical parameters: uniform spacing h, diffusivity D (heat equation du/dt = D * Laplacian(u))
 * - Optional `stored_dt` (demo / documentation only; not used by `DiffusionSolver`)
 */
struct Diffusion3DContext
{
    int nx, ny, nz; // Grid dimensions (x, y, z)

    /**
     * @brief GPU block dimensions (threads per block)
     */
    int bx = 8;
    int by = 8;
    int bz = 8;

    int n; // Total number of grid points (nx * ny * nz)

    std::vector<double> u;      // Current scalar field (state at time t)
    std::vector<double> u_next; // Scratch buffer for explicit Euler (same layout as `u`)

    /**
     * @brief Optional reference timestep (e.g. macro tick for demos). Not read by `DiffusionSolver`.
     */
    double stored_dt;

    /**
     * @brief Uniform grid spacing (same units as the physical domain).
     */
    double h;

    /**
     * @brief Diffusivity D in du/dt = D * Laplacian(u); units [length]^2/[time] with consistent h, dt.
     */
    double D;

private:
    enum class Init
    {
        kInternal
    };

    Diffusion3DContext(Init,
                       int nx,
                       int ny,
                       int nz,
                       double h,
                       double D,
                       double stored_dt)
        : nx(nx),
          ny(ny),
          nz(nz),
          n(nx * ny * nz),
          u(static_cast<size_t>(n), 0.0),
          u_next(static_cast<size_t>(n), 0.0),
          stored_dt(stored_dt),
          h(h),
          D(D)
    {
        assert(h > 0.0);
        assert(D >= 0.0);
    }

public:
    /**
     * @brief Grid + physics only; `stored_dt` is set to 0 (unused). Preferred for tests and library use.
     */
    static Diffusion3DContext make(int nx, int ny, int nz, double h, double D)
    {
        return Diffusion3DContext(Init::kInternal, nx, ny, nz, h, D, 0.0);
    }

    /**
     * @brief Same as `make`, plus a stored reference timestep (e.g. `tick_dt` in demos). Not used by `DiffusionSolver`.
     */
    static Diffusion3DContext make_with_reference_dt(int nx, int ny, int nz, double h, double D, double reference_dt)
    {
        return Diffusion3DContext(Init::kInternal, nx, ny, nz, h, D, reference_dt);
    }

    /**
     * @brief Convert 3D coordinates (x,y,z) into a flattened 1D index.
     *
     * index = (z * ny + y) * nx + x
     */
    inline int idx(int x, int y, int z) const
    {
        return (z * ny + y) * nx + x;
    }

    /**
     * Modify field value at a coordinate.
     */
    void set_at_coord(int x, int y, int z, double value)
    {
        u[idx(x, y, z)] = value;
    }
};

/**
 * @brief One explicit Euler step on the 6-point stencil with micro-step `dt_sub`.
 *
 * Reads `ctx.u`, writes `ctx.u_next`, then swaps buffers so the updated state lives in `ctx.u`.
 * Neighbors outside the data domain are omitted.
 *
 * @param dt_sub Time increment for this step (often `tick_dt / n_sub`).
 */
void diffusion3d_step_euler_cpu(Diffusion3DContext &ctx, double dt_sub);

/**
 * @brief One GPU explicit Euler step on the 6-point stencil (device buffers).
 *
 * Launches `diffusion3d_kernel` with micro-step `dt_sub`, synchronizes, then swaps `d_u` and
 * `d_u_next` so the updated field is referenced by `d_u`.
 *
 * @param dt_sub Same role as for `diffusion3d_step_euler_cpu`.
 * @param d_u,d_u_next Device pointers to fields of length `ctx.n`; swapped after the step.
 */
void diffusion3d_step_euler_gpu(Diffusion3DContext &ctx, double dt_sub, double *&d_u, double *&d_u_next);

#endif // DIFFUSION3D_COMMON_H
