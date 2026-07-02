#ifndef DIFFUSION3D_STEP_EULER_CPU_H
#define DIFFUSION3D_STEP_EULER_CPU_H

/**
 * @file diffusion3d_step_euler_cpu.h
 * @brief CPU explicit Euler plus 6-point face Laplacian (open box).
 */

#include "scalar_field_grid.h"

/**
 * @brief One explicit Euler micro-step on `grid` (open-box 6-point stencil).
 */
void diffusion3d_step_euler_scalar(ScalarFieldGrid &grid, double h, double D, double dt_sub);

#endif // DIFFUSION3D_STEP_EULER_CPU_H
