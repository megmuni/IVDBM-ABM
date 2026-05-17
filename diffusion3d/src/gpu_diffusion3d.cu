#include <cassert>
#include <utility>

#include "gpu_diffusion3d.cuh"

#ifdef __CUDACC__
__device__ inline bool is_inside_domain(int x, int y, int z, int nx, int ny, int nz)
{
    return x >= 0 && x < nx &&
           y >= 0 && y < ny &&
           z >= 0 && z < nz;
}
/**
 * @brief GPU kernel for one timestep of 3D diffusion (explicit Euler + 6-point stencil, no convolution).
 *
 * Each CUDA thread updates one patch using its immediate neighbors.
 */
__global__ void diffusion3d_kernel(
    double *d_u,
    double *d_u_next,
    int nx, int ny, int nz,
    double lap_scale,
    double dt)
{
    /**
     * Extract 3D cordinates from block and thread indices.
     * Each thread in each block corresponds to one patch, with each
     * thread computing that patch's new value based on its current value.
     */
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x >= nx || y >= ny || z >= nz) return;

    int idx = (z * ny + y) * nx + x; // Flattened index for 3D coordinates (x,y,z)

    double center = d_u[idx]; // Current value at the patch being updated
    double lap = 0.0; // Accumulator for Laplacian contribution from neighbors

    /**
     * Direction offsets for 6-point stencil (x±1, y±1, z±1).
     * {1, 0, 0} means +1 in x, {0, 1, 0} means +1 in y, etc.
     * So:
     * - { 1, 0, 0} --> right neighbor
     * - {-1, 0, 0} --> left neighbor
     * - { 0, 1, 0} --> forward neighbor
     * - { 0,-1, 0} --> backward neighbor
     * - { 0, 0, 1} --> up neighbor
     * - { 0, 0,-1} --> down neighbor
     */
    const int nbr[6][3] = {
        { 1, 0, 0}, {-1, 0, 0},
        { 0, 1, 0}, { 0,-1, 0},
        { 0, 0, 1}, { 0, 0,-1}
    };

    /**
     * Iterate over 6 neighbours, compute contribution to
     * Laplacian (neighbor value - center value) and accumulate.
     */
    for (int k = 0; k < 6; ++k)
    {
        int nx_ = x + nbr[k][0];
        int ny_ = y + nbr[k][1];
        int nz_ = z + nbr[k][2];

        if (is_inside_domain(nx_, ny_, nz_, nx, ny, nz))
        {
            int nidx = (nz_ * ny + ny_) * nx + nx_;
            lap += d_u[nidx] - center; // Contribution to Laplacian
        }
    }

    /**
     * Update rule (same as CPU): u_next = u + dt * (D/h^2) * sum_face (u_nei - u).
     * lap_scale is D/h^2 computed on host.
     */
    d_u_next[idx] = center + lap_scale * dt * lap;
}

#endif

/**
 * @brief Host wrapper: launch `diffusion3d_kernel`, sync, swap device buffers so `d_u` holds the new field.
 */
void diffusion3d_step_euler_gpu(Diffusion3DContext &ctx, double dt_sub, double *&d_u, double *&d_u_next)
{
#ifdef __CUDACC__
    assert(ctx.h > 0.0);
    assert(ctx.D >= 0.0);
    assert(dt_sub >= 0.0);

    const double lap_scale = ctx.D / (ctx.h * ctx.h);

    dim3 block(ctx.bx, ctx.by, ctx.bz);
    dim3 grid(
        (ctx.nx + ctx.bx - 1) / ctx.bx,
        (ctx.ny + ctx.by - 1) / ctx.by,
        (ctx.nz + ctx.bz - 1) / ctx.bz
    );

    diffusion3d_kernel<<<grid, block>>>(
        d_u,
        d_u_next,
        ctx.nx,
        ctx.ny,
        ctx.nz,
        lap_scale,
        dt_sub
    );

    cudaDeviceSynchronize();
    std::swap(d_u, d_u_next);
#else
    (void)ctx;
    (void)dt_sub;
    (void)d_u;
    (void)d_u_next;
#endif
}