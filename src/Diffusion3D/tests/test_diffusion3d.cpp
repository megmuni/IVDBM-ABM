/**
 * @file test_diffusion3d.cpp
 * @brief Catch2 tests for Diffusion3D Phase I core.
 */

#include <catch2/catch.hpp>

#include <cmath>
#include <vector>

#include "diffusion3d_timestep.h"
#include "diffusion3d_step_euler_cpu.h"
#include "explicit_multi_species_heat_stepper.h"
#include "explicit_gpu_stencil_stepper.h"
#include "make_multi_species_diffusion_engine.h"
#include "multi_species_diffusion_settings.h"
#include "multi_species_field_grid.h"

namespace
{

bool nearly_equal_vec(const std::vector<double> &a, const std::vector<double> &b, double eps)
{
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
    {
        if (std::fabs(a[i] - b[i]) > eps)
            return false;
    }
    return true;
}

void fill_ramp(MultiSpeciesFieldGrid &grid, SpeciesId id)
{
    ScalarFieldGrid &g = grid.grid(id);
    for (int z = 0; z < g.nz_; ++z)
        for (int y = 0; y < g.ny_; ++y)
            for (int x = 0; x < g.nx_; ++x)
                g.at(x, y, z) = 0.01 * static_cast<double>(x + y * g.nx_ + z * g.nx_ * g.ny_);
}

MultiSpeciesDiffusionSettings make_settings(SpeciesId id, double D, double safety = 1.0)
{
    MultiSpeciesDiffusionSettings settings;
    settings.species_diffusivities[id] = D;
    settings.safety = safety;
    settings.algorithm = DiffusionAlgorithm::ExplicitHeatEquation;
    return settings;
}

} // namespace

TEST_CASE("MultiSpeciesFieldGrid holds per-species grids on shared domain", "[core][field_grid]")
{
    const std::vector<SpeciesId> ids = {0, 1, 2};
    MultiSpeciesFieldGrid grid(ids, 4, 3, 2, 0.5);

    CHECK(grid.nx_ == 4);
    CHECK(grid.ny_ == 3);
    CHECK(grid.nz_ == 2);
    CHECK(grid.h_ == Approx(0.5));
    CHECK(grid.species().size() == 3);
    CHECK(grid.has_species(1));

    grid.grid(0).at(1, 2, 1) = 3.5;
    CHECK(grid.grid(0).at(1, 2, 1) == Approx(3.5));
    CHECK(grid.grid(1).at(1, 2, 1) == Approx(0.0));
}

TEST_CASE("Single-species stepper matches explicit scalar euler micro-steps", "[core][engine]")
{
    const int nx = 5, ny = 4, nz = 3;
    const double h = 1.0;
    const double D = 1.0;
    const double tick_dt = 0.5;
    const SpeciesId species = 0;

    const SubstepPlan plan = plan_substeps(tick_dt, compute_stability_constraint(h, D, 1.0));

    ScalarFieldGrid ref(nx, ny, nz);
    MultiSpeciesFieldGrid grid({species}, nx, ny, nz, h);
    fill_ramp(grid, species);
    ref.data_ = grid.grid(species).data_;

    for (int i = 0; i < plan.n_sub; ++i)
        diffusion3d_step_euler_scalar(ref, h, D, plan.dt_sub);

    ExplicitMultiSpeciesHeatStepper stepper;
    MultiSpeciesDiffusionSettings settings = make_settings(species, D);
    stepper.configure_species_interval(grid, settings, tick_dt);
    stepper.advance_species_interval(grid);

    CHECK(nearly_equal_vec(ref.data_, grid.grid(species).data_, 1e-13));
}

TEST_CASE("Multi-species fields evolve independently", "[core][engine][multi]")
{
    const int nx = 4, ny = 4, nz = 4;
    const double h = 1.0;
    const double tick_dt = 0.2;
    const SpeciesId slow = 0;
    const SpeciesId fast = 1;

    MultiSpeciesFieldGrid grid({slow, fast}, nx, ny, nz, h);
    fill_ramp(grid, slow);
    fill_ramp(grid, fast);

    const std::vector<double> initial_slow = grid.grid(slow).data_;
    const std::vector<double> initial_fast = grid.grid(fast).data_;

    MultiSpeciesDiffusionSettings settings;
    settings.species_diffusivities[slow] = 0.1;
    settings.species_diffusivities[fast] = 1.0;
    settings.safety = 1.0;
    settings.algorithm = DiffusionAlgorithm::ExplicitHeatEquation;

    ExplicitMultiSpeciesHeatStepper stepper;
    stepper.configure_species_interval(grid, settings, tick_dt);
    REQUIRE(stepper.n_sub(fast) >= stepper.n_sub(slow));
    stepper.advance_species_interval(grid);

    CHECK_FALSE(nearly_equal_vec(initial_slow, grid.grid(slow).data_, 1e-15));
    CHECK_FALSE(nearly_equal_vec(initial_fast, grid.grid(fast).data_, 1e-15));
    CHECK_FALSE(nearly_equal_vec(grid.grid(slow).data_, grid.grid(fast).data_, 1e-12));
}

TEST_CASE("MakeMultiSpeciesDiffusionEngine creates CPU stepper", "[core][factory]")
{
    MultiSpeciesDiffusionSettings settings = make_settings(0, 0.5);
    auto engine = MakeMultiSpeciesDiffusionEngine(settings);
    REQUIRE(engine != nullptr);
    CHECK(engine->resolved_algorithm() == DiffusionAlgorithm::ExplicitHeatEquation);
}

TEST_CASE("configure_species_interval requires tick_dt macro interval", "[core][engine]")
{
    const SpeciesId id = 0;
    MultiSpeciesFieldGrid grid({id}, 3, 3, 3, 1.0);
    fill_ramp(grid, id);
    const std::vector<double> initial = grid.grid(id).data_;

    ExplicitMultiSpeciesHeatStepper stepper;
    MultiSpeciesDiffusionSettings settings = make_settings(id, 1.0);

    stepper.configure_species_interval(grid, settings, 0.0);
    stepper.advance_species_interval(grid);
    CHECK(nearly_equal_vec(initial, grid.grid(id).data_, 1e-15));

    MultiSpeciesFieldGrid grid2({id}, 3, 3, 3, 1.0);
    fill_ramp(grid2, id);
    const auto before = grid2.grid(id).data_;

    stepper.configure_species_interval(grid2, settings, 0.5);
    stepper.advance_species_interval(grid2);
    CHECK_FALSE(nearly_equal_vec(before, grid2.grid(id).data_, 1e-15));
}

#ifdef DIFFUSION3D_CUDA

TEST_CASE("GpuStencil stepper matches CPU stepper for single species", "[core][gpu]")
{
    const int nx = 6, ny = 5, nz = 4;
    const double h = 0.4, D = 0.08, tick_dt = 0.07;
    const SpeciesId species = 0;

    MultiSpeciesFieldGrid cpu_grid({species}, nx, ny, nz, h);
    MultiSpeciesFieldGrid gpu_grid({species}, nx, ny, nz, h);
    fill_ramp(cpu_grid, species);
    fill_ramp(gpu_grid, species);

    MultiSpeciesDiffusionSettings settings = make_settings(species, D);

    ExplicitMultiSpeciesHeatStepper cpu_stepper;
    ExplicitGpuStencilStepper gpu_stepper;
    cpu_stepper.configure_species_interval(cpu_grid, settings, tick_dt);
    gpu_stepper.configure_species_interval(gpu_grid, settings, tick_dt);
    cpu_stepper.advance_species_interval(cpu_grid);
    gpu_stepper.advance_species_interval(gpu_grid);

    CHECK(nearly_equal_vec(cpu_grid.grid(species).data_, gpu_grid.grid(species).data_, 1e-12));
}

TEST_CASE("MakeMultiSpeciesDiffusionEngine creates GpuStencil stepper", "[core][factory][gpu]")
{
    MultiSpeciesDiffusionSettings settings = make_settings(0, 0.5);
    settings.algorithm = DiffusionAlgorithm::GpuStencil;
    auto engine = MakeMultiSpeciesDiffusionEngine(settings);
    REQUIRE(engine != nullptr);
    CHECK(engine->resolved_algorithm() == DiffusionAlgorithm::GpuStencil);
}

#endif // DIFFUSION3D_CUDA
