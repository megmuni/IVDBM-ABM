/**
 * @file diffusion3d_timestep_test.cpp
 * @brief Catch2 tests for diffusion3d_timestep (stability constraint + plan_substeps).
 */

#include <cmath>
#include <limits>

#include <catch2/catch.hpp>

#include "diffusion3d_timestep.h"

TEST_CASE("compute_stability_constraint matches h^2/(6D) with safety", "[timestep]")
{
    const double h = 0.5;
    const double D = 0.12;
    const double safety = 1.0;
    const double expected = safety * h * h / (6.0 * D);
    CHECK(compute_stability_constraint(h, D, safety) == Approx(expected));
}

TEST_CASE("compute_stability_constraint scales with safety", "[timestep]")
{
    const double h = 1.0;
    const double D = 0.1;
    const double full = compute_stability_constraint(h, D, 1.0);
    const double half = compute_stability_constraint(h, D, 0.5);
    CHECK(half == Approx(0.5 * full));
}

TEST_CASE("compute_stability_constraint is infinite when D is zero", "[timestep]")
{
    const double dt_max = compute_stability_constraint(1.0, 0.0, 1.0);
    CHECK(std::isinf(dt_max));
    CHECK(dt_max > 0.0);
}

TEST_CASE("plan_substeps tick_dt zero yields one zero-length substep", "[timestep]")
{
    const SubstepPlan p = plan_substeps(0.0, 1e-3);
    CHECK(p.n_sub == 1);
    CHECK(p.dt_sub == Approx(0.0));
}

TEST_CASE("plan_substeps with infinite dt_max uses single step", "[timestep]")
{
    const SubstepPlan p = plan_substeps(2.5, std::numeric_limits<double>::infinity());
    CHECK(p.n_sub == 1);
    CHECK(p.dt_sub == Approx(2.5));
}

TEST_CASE("plan_substeps ceil splits tick into substeps bounded by dt_max", "[timestep]")
{
    const double tick_dt = 1.0;
    const double dt_max = 0.3;
    const SubstepPlan p = plan_substeps(tick_dt, dt_max);
    CHECK(p.n_sub == 4);
    CHECK(p.dt_sub == Approx(0.25));
    CHECK(p.dt_sub <= dt_max * (1.0 + 1e-12));
}

TEST_CASE("plan_substeps exact division yields integer ratio", "[timestep]")
{
    const SubstepPlan p = plan_substeps(1.0, 0.25);
    CHECK(p.n_sub == 4);
    CHECK(p.dt_sub == Approx(0.25));
}
