#ifndef DIFFUSION3D_GPU_H
#define DIFFUSION3D_GPU_H

/**
 * @file diffusion3d_gpu.h
 * @brief GPU explicit Euler stencil kernel (device buffers).
 */

/**
 * @brief One GPU explicit Euler step on the 6-point stencil.
 *
 * Launches the diffusion kernel, synchronizes, then swaps `d_u` and `d_u_next`
 * so the updated field is referenced by `d_u`.
 *
 * @param nx, ny, nz Grid dimensions
 * @param h Uniform spacing, @param D diffusivity
 * @param dt_sub Explicit Euler micro-step
 * @param d_u, d_u_next Device buffers of length nx*ny*nz (swapped after step)
 * @param bx, by, bz CUDA block dimensions (default 8)
 */
void diffusion3d_step_euler_gpu(int nx, int ny, int nz,
                                double h, double D, double dt_sub,
                                double *&d_u, double *&d_u_next,
                                int bx = 8, int by = 8, int bz = 8);

#endif // DIFFUSION3D_GPU_H
