#pragma once

#include <vector>

/**
 * @file analytic_diffusion3d.h
 * @brief Analytic reference for 3D heat/diffusion from a point impulse (test support).
 *
 *   du/dt = D * Laplacian(u)
 *
 * Row-major layout: idx = x + nx*(y + ny*z) (same as ScalarFieldGrid).
 */

std::vector<double> analytic_gaussian_impulse_3d(int nx, int ny, int nz,
                                                   double h, double D,
                                                   double t,
                                                   int cx, int cy, int cz,
                                                   double mass = 1.0);

double rel_l2_error(const std::vector<double> &a, const std::vector<double> &b);
