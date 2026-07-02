/**
 * @file engine_stepper_test.cpp
 * @brief Catch2 tests for MultiSpeciesDiffusionEngine steppers (replaces DiffusionSolver suite).
 */

#include <catch2/catch.hpp>

#include <cmath>
#include <vector>

#include "diffusion3d_step_euler_cpu.h"
#include "diffusion3d_timestep.h"
#include "explicit_gpu_fft_stepper.h"
#include "explicit_multi_species_heat_stepper.h"
#include "make_multi_species_diffusion_engine.h"
#include "multi_species_diffusion_settings.h"
#include "multi_species_field_grid.h"

namespace
{

bool nearly_equal_vec(const std::vector<double> &a, const std::vector<double> &b, double eps)
{
    REQUIRE(a.size() == b.size());
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

MultiSpeciesDiffusionSettings make_settings(SpeciesId id, double D, DiffusionAlgorithm algo,
                                            double safety = 1.0)
{
    MultiSpeciesDiffusionSettings settings;
    settings.species_diffusivities[id] = D;
    settings.safety = safety;
    settings.algorithm = algo;
    return settings;
}

} // namespace

TEST_CASE("ExplicitMultiSpeciesHeatStepper splits tick into stable substeps", "[engine][cpu]")
{
    const int nx = 5, ny = 4, nz = 3;
    const double h = 1.0;
    const double D = 1.0;
    const double safety = 1.0;
    const double tick_dt = 0.5;
    const SpeciesId species = 0;

    const double dt_max = compute_stability_constraint(h, D, safety);
    const SubstepPlan plan = plan_substeps(tick_dt, dt_max);
    REQUIRE(plan.n_sub >= 2);

    ScalarFieldGrid ref(nx, ny, nz);
    MultiSpeciesFieldGrid grid({species}, nx, ny, nz, h);
    fill_ramp(grid, species);
    ref.data_ = grid.grid(species).data_;

    for (int i = 0; i < plan.n_sub; ++i)
        diffusion3d_step_euler_scalar(ref, h, D, plan.dt_sub);

    ExplicitMultiSpeciesHeatStepper stepper;
    MultiSpeciesDiffusionSettings settings = make_settings(species, D, DiffusionAlgorithm::ExplicitHeatEquation, safety);
    stepper.configure_species_interval(grid, settings, tick_dt);
    stepper.advance_species_interval(grid);

    CHECK(nearly_equal_vec(ref.data_, grid.grid(species).data_, 1e-13));
}

TEST_CASE("MakeMultiSpeciesDiffusionEngine resolves algorithm per CUDA build", "[engine][backend]")
{
    MultiSpeciesDiffusionSettings settings = make_settings(0, 0.5, DiffusionAlgorithm::GpuFftPrecomputed);
    auto engine = MakeMultiSpeciesDiffusionEngine(settings);
    REQUIRE(engine != nullptr);
#ifndef DIFFUSION3D_CUDA
    CHECK(engine->resolved_algorithm() == DiffusionAlgorithm::ExplicitHeatEquation);
#else
    CHECK(engine->resolved_algorithm() == DiffusionAlgorithm::GpuFftPrecomputed);
#endif
}

#ifdef DIFFUSION3D_CUDA

TEST_CASE("ExplicitGpuFftStepper matches CPU stepper in interior (small tick)", "[engine][gpu][fft]")
{
    const int nx = 6, ny = 6, nz = 6;
    const double tick_dt = 1e-3;
    const double h = 1.0;
    const double D = 0.1;
    const SpeciesId species = 0;

    MultiSpeciesFieldGrid cpu_grid({species}, nx, ny, nz, h);
    MultiSpeciesFieldGrid fft_grid({species}, nx, ny, nz, h);

    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const double v =
                    0.01 * static_cast<double>((z + 1) * 100 + (y + 1) * 10 + (x + 1));
                cpu_grid.grid(species).at(x, y, z) = v;
                fft_grid.grid(species).at(x, y, z) = v;
            }
        }
    }

    MultiSpeciesDiffusionSettings cpu_settings = make_settings(species, D, DiffusionAlgorithm::ExplicitHeatEquation);
    MultiSpeciesDiffusionSettings fft_settings = make_settings(species, D, DiffusionAlgorithm::GpuFftPrecomputed);
    fft_settings.fft_real_extent_x = 8;
    fft_settings.fft_real_extent_y = 8;
    fft_settings.fft_real_extent_z = 8;

    ExplicitMultiSpeciesHeatStepper cpu_stepper;
    ExplicitGpuFftStepper fft_stepper;
    cpu_stepper.configure_species_interval(cpu_grid, cpu_settings, tick_dt);
    fft_stepper.configure_species_interval(fft_grid, fft_settings, tick_dt);
    cpu_stepper.advance_species_interval(cpu_grid);
    fft_stepper.advance_species_interval(fft_grid);

    const ScalarFieldGrid &cpu = cpu_grid.grid(species);
    const ScalarFieldGrid &fft = fft_grid.grid(species);

    double max_abs = 0.0;
    double num = 0.0;
    double den = 0.0;
    for (int z = 1; z < nz - 1; ++z)
    {
        for (int y = 1; y < ny - 1; ++y)
        {
            for (int x = 1; x < nx - 1; ++x)
            {
                const double cu = cpu.at(x, y, z);
                const double fu = fft.at(x, y, z);
                const double diff = fu - cu;
                max_abs = std::max(max_abs, std::fabs(diff));
                num += diff * diff;
                den += cu * cu;
            }
        }
    }
    const double rel_l2 = (den > 0.0) ? std::sqrt(num / den) : std::sqrt(num);
    INFO("rel_l2=" << rel_l2 << " max_abs=" << max_abs);
    CHECK(rel_l2 < 1e-4);
    CHECK(max_abs < 1e-3);
}

#endif // DIFFUSION3D_CUDA
