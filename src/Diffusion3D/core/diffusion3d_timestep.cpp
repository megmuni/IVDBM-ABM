/**
 * @file diffusion3d_timestep.cpp
 * @brief Implements stability bound and substep planning for explicit Euler + 6-point 3D Laplacian.
 */

#include "diffusion3d_timestep.h"

#include <cassert>
#include <cmath>
#include <limits>

double compute_stability_constraint(double h, double D, double safety)
{
    assert(h > 0.0);
    assert(D >= 0.0);
    assert(safety > 0.0 && safety <= 1.0);

    if (D == 0.0)
        return std::numeric_limits<double>::infinity();

    return safety * (h * h) / (6.0 * D);
}

SubstepPlan plan_substeps(double tick_dt, double dt_max)
{
    assert(tick_dt >= 0.0);

    if (tick_dt == 0.0)
        return SubstepPlan{1, 0.0};

    if (!std::isfinite(dt_max))
        return SubstepPlan{1, tick_dt};

    assert(dt_max > 0.0);

    const double ratio = tick_dt / dt_max;
    const int n_sub = std::max(1, static_cast<int>(std::ceil(ratio)));
    const double dt_sub = tick_dt / static_cast<double>(n_sub);
    return SubstepPlan{n_sub, dt_sub};
}
