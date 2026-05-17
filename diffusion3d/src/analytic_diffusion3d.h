#pragma once

#include "diffusion3d_common.h"
#include <vector>

/**
 * @file analytic_diffusion3d.h
 * @brief Analytic reference for 3D heat/diffusion from a point impulse.
 *
 * Provides a simple Gaussian Green's-function reference for:
 *
 *   du/dt = D * Laplacian(u)
 *
 * with an impulse initial condition. Grid indices are converted to physical
 * separation using ctx.h (cell spacing).
 */

/**
 * @brief Compute the analytic solution u(x,y,z,t) for an impulse at (cx,cy,cz).
 *
 * Discrete sampling of the continuous Green's function:
 *
 *   u(r,t) = M / (4π D t)^(3/2) * exp(-r^2 / (4 D t))
 *
 * where M is the total impulse mass.
 */
std::vector<double> analytic_gaussian_impulse_3d(const Diffusion3DContext &ctx,
                                                 double t,
                                                 int cx,
                                                 int cy,
                                                 int cz,
                                                 double mass = 1.0);

/**
 * @brief Compute relative L2 error between two equally-sized fields.
 */
double rel_l2_error(const std::vector<double> &a, const std::vector<double> &b);

