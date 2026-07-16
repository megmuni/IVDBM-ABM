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

TEST_CASE("stencil centre weight stays positive for every D at safety < 1", "[timestep]")
{
    // The 6-point update weights the centre cell by 1 - 6s, s = D*dt_sub/h^2.
    // safety < 1 must keep that weight >= 1 - safety > 0 regardless of D, h or
    // tick_dt; at safety == 1 it collapses to 0 and the cell is replaced by the
    // mean of its neighbours (even/odd sublattices decouple -> checkerboard).
    //
    // The biological D/h/tick below are the dangerous case: 6*D*tick/h^2 is an
    // exact integer (324 for TNF, 39096 for O2), so plan_substeps divides evenly
    // and dt_sub lands exactly on dt_max rather than safely under it.
    const double h = 0.01;
    const double tick_dt = 30.0;
    const double safety = 0.5;

    for (const double D : {0.00018, 0.02172, 0.1, 1.0, 1e-6})
    {
        const double dt_max = compute_stability_constraint(h, D, safety);
        const SubstepPlan plan = plan_substeps(tick_dt, dt_max);
        const double s = D * plan.dt_sub / (h * h);
        const double centre_weight = 1.0 - 6.0 * s;

        INFO("D=" << D << " n_sub=" << plan.n_sub << " s=" << s
                  << " centre_weight=" << centre_weight);
        CHECK(s <= safety / 6.0 + 1e-12);
        CHECK(centre_weight >= 1.0 - safety - 1e-12);
        CHECK(centre_weight > 0.0);
    }
}
