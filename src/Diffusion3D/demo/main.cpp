/**
 * @file main.cpp
 * @brief ParaView demo: CPU stencil, GPU stencil, and GPU FFT via
 * MultiSpeciesDiffusionEngine.
 *
 * Writes time-series `.vti` under output/cpu_series/, output/gpu_series/, and
 * output/fft_series/ (CUDA build). Creates output/ if needed.
 * Open in ParaView: File → Open → select diffusion3d_t*.vti → group as time
 * series.
 */

#include <cassert>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>

#include "diffusion3d_timestep.h"
#include "explicit_multi_species_heat_stepper.h"
#include "make_multi_species_diffusion_engine.h"
#include "multi_species_diffusion_settings.h"
#include "multi_species_field_grid.h"
#include "vtk_export.h"

namespace {

constexpr SpeciesId kSpecies = 0;
constexpr const char *kOutputRoot = "output";

const char *algorithm_label(DiffusionAlgorithm algo) {
  switch (algo) {
  case DiffusionAlgorithm::ExplicitHeatEquation:
    return "ExplicitHeatEquation (CPU stencil)";
  case DiffusionAlgorithm::GpuStencil:
    return "GpuStencil";
  case DiffusionAlgorithm::GpuFftPrecomputed:
    return "GpuFftPrecomputed (FFT)";
  }
  return "unknown";
}

MultiSpeciesFieldGrid make_initial_grid(int nx, int ny, int nz, double h) {
  MultiSpeciesFieldGrid grid({kSpecies}, nx, ny, nz, h);
  grid.grid(kSpecies).at(nx / 2, ny / 2, nz / 2) = 1.0;
  return grid;
}

MultiSpeciesDiffusionSettings make_settings(DiffusionAlgorithm algo, double D) {
  MultiSpeciesDiffusionSettings settings;
  settings.species_diffusivities[kSpecies] = D;
  settings.safety = 1.0;
  settings.algorithm = algo;
  if (algo == DiffusionAlgorithm::GpuFftPrecomputed) {
    settings.fft_real_extent_x = 32;
    settings.fft_real_extent_y = 32;
    settings.fft_real_extent_z = 32;
  }
  return settings;
}

void log_substep_plan(const char *label, double h, double D, double tick_dt,
                      double safety) {
  const double dt_max = compute_stability_constraint(h, D, safety);
  const SubstepPlan plan = plan_substeps(tick_dt, dt_max);
  std::cout << "demo [" << label << "] plan: n_sub=" << plan.n_sub
            << " dt_sub=" << plan.dt_sub << " (tick_dt=" << tick_dt << ")\n";
}

void run_time_series(const char *label, const std::string &outdir,
                     DiffusionAlgorithm algo, int nx, int ny, int nz, double h,
                     double D, double tick_dt, int nsteps) {
  std::error_code ec;
  std::filesystem::create_directories(outdir, ec);
  assert(!ec);

  MultiSpeciesFieldGrid grid = make_initial_grid(nx, ny, nz, h);
  MultiSpeciesDiffusionSettings settings = make_settings(algo, D);

  auto engine = MakeMultiSpeciesDiffusionEngine(settings);
  engine->configure_species_interval(grid, settings, tick_dt);

  std::cout << "demo [" << label << "] "
            << algorithm_label(engine->resolved_algorithm()) << "\n";
  log_substep_plan(label, h, D, tick_dt, settings.safety);

  if (auto *cpu =
          dynamic_cast<ExplicitMultiSpeciesHeatStepper *>(engine.get())) {
    assert(cpu->n_sub(kSpecies) > 1 &&
           "pick tick_dt/h/D that triggers substepping");
  }

  for (int t = 0; t <= nsteps; ++t) {
    char path[256];
    std::snprintf(path, sizeof(path), "%s/diffusion3d_t%d.vti", outdir.c_str(),
                  t);
    export_scalar_field_to_vti(grid.grid(kSpecies), h, path, "u");
    if (t < nsteps)
      engine->advance_species_interval(grid);
  }

  std::cout << "demo [" << label << "] wrote " << (nsteps + 1) << " frames to "
            << outdir << "/\n";
}

#ifdef DIFFUSION3D_CUDA
void export_cpu_gpu_stencil_compare(int nx, int ny, int nz, double h, double D,
                                    double tick_dt, int nsteps) {
  const std::string outdir = std::string(kOutputRoot) + "/compare_series";
  std::error_code ec;
  std::filesystem::create_directories(outdir, ec);
  assert(!ec);

  MultiSpeciesFieldGrid cpu_grid = make_initial_grid(nx, ny, nz, h);
  MultiSpeciesFieldGrid gpu_grid = make_initial_grid(nx, ny, nz, h);

  MultiSpeciesDiffusionSettings cpu_settings =
      make_settings(DiffusionAlgorithm::ExplicitHeatEquation, D);
  MultiSpeciesDiffusionSettings gpu_settings =
      make_settings(DiffusionAlgorithm::GpuStencil, D);

  auto cpu_engine = MakeMultiSpeciesDiffusionEngine(cpu_settings);
  auto gpu_engine = MakeMultiSpeciesDiffusionEngine(gpu_settings);
  cpu_engine->configure_species_interval(cpu_grid, cpu_settings, tick_dt);
  gpu_engine->configure_species_interval(gpu_grid, gpu_settings, tick_dt);

  for (int t = 0; t < nsteps; ++t) {
    cpu_engine->advance_species_interval(cpu_grid);
    gpu_engine->advance_species_interval(gpu_grid);
  }

  export_cpu_gpu_compare_to_vti(cpu_grid.grid(kSpecies).data_,
                                gpu_grid.grid(kSpecies).data_, nx, ny, nz, h,
                                outdir + "/cpu_vs_gpu_stencil_final.vti");

  std::cout << "demo [compare] wrote " << outdir
            << "/cpu_vs_gpu_stencil_final.vti\n";
}
#endif

} // namespace

int main() {
  const int nsteps = 4;
  const int nx = 8, ny = 8, nz = 8;
  const double tick_dt = 2.0;
  const double h = 1.0;
  const double D = 0.1;

  std::cout << "diffusion3d_demo: nx=" << nx << " ny=" << ny << " nz=" << nz
            << " h=" << h << " D=" << D << " tick_dt=" << tick_dt
            << " nsteps=" << nsteps << "\n";

  std::error_code ec;
  std::filesystem::create_directories(kOutputRoot, ec);
  assert(!ec);

  run_time_series("cpu", std::string(kOutputRoot) + "/cpu_series",
                  DiffusionAlgorithm::ExplicitHeatEquation, nx, ny, nz, h, D,
                  tick_dt, nsteps);

#ifdef DIFFUSION3D_CUDA
  run_time_series("gpu_stencil", std::string(kOutputRoot) + "/gpu_series",
                  DiffusionAlgorithm::GpuStencil, nx, ny, nz, h, D, tick_dt,
                  nsteps);

  run_time_series("gpu_fft", std::string(kOutputRoot) + "/fft_series",
                  DiffusionAlgorithm::GpuFftPrecomputed, nx, ny, nz, h, D,
                  tick_dt, nsteps);

  export_cpu_gpu_stencil_compare(nx, ny, nz, h, D, tick_dt, nsteps);
#else
  std::cout
      << "demo: GPU paths skipped - rebuild with -DDIFFUSION3D_CUDA=ON for "
         "output/gpu_series/, output/fft_series/, and output/compare_series/\n";
#endif

  std::cout << "diffusion3d_demo: done. Open *.vti in ParaView.\n";
  return 0;
}
