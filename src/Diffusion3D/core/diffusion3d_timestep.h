#ifndef DIFFUSION3D_TIMESTEP_H
#define DIFFUSION3D_TIMESTEP_H

/**
 * @file diffusion3d_timestep.h
 * @brief CFL-style bound (von Neumann, h^2/(6D)) and macro-to-micro timestep splitting.
 */

/**
 * @brief Result of splitting a macro interval `tick_dt` into stable explicit Euler micro-steps.
 */
struct SubstepPlan
{
    int n_sub;     ///< Number of micro-steps (>= 1).
    double dt_sub; ///< tick_dt / n_sub (non-negative).
};

/**
 * @brief Maximum admissible explicit Euler step for the 3D 6-point stencil (diffusion only).
 *
 * @param h Uniform grid spacing (> 0).
 * @param D Diffusivity (>= 0). If D == 0, returns +infinity (no diffusion restriction).
 * @param safety Factor in (0, 1] applied to the theoretical bound h^2/(6D).
 */
double compute_stability_constraint(double h, double D, double safety);

/**
 * @brief Choose n_sub and dt_sub so n_sub * dt_sub == tick_dt and each micro-step stays within dt_max
 * when dt_max is finite and positive.
 *
 * tick_dt == 0: returns { 1, 0 }.
 * dt_max non-finite (NaN, +/-infinity): returns { 1, tick_dt } (no splitting).
 */
SubstepPlan plan_substeps(double tick_dt, double dt_max);

#endif
