/**
 * @file explicit_gpu_fft_stepper.cpp
 * @brief GPU FFT multi-species stepper.
 */

#include "explicit_gpu_fft_stepper.h"
#include "diffusion3d_fft_scratch.h"
#include <stdexcept>

ExplicitGpuFftStepper::ExplicitGpuFftStepper() = default;

ExplicitGpuFftStepper::~ExplicitGpuFftStepper() {
  for (auto &[id, scratch] : fft_scratch_) {
    (void)id;
    diffusion3d_fft_scratch_destroy(scratch);
  }
  fft_scratch_.clear();
}

void ExplicitGpuFftStepper::destroy_fft_scratch(SpeciesId id) {
  auto it = fft_scratch_.find(id);
  if (it != fft_scratch_.end()) {
    diffusion3d_fft_scratch_destroy(it->second);
    fft_scratch_.erase(it);
  }
}

static int derive_fft_axis_len(int domain_len) { return 2 * domain_len; }

void ExplicitGpuFftStepper::configure_species_interval(
    const MultiSpeciesFieldGrid &grid,
    const MultiSpeciesDiffusionSettings &settings, double tick_dt) {
#ifndef DIFFUSION3D_CUDA
  throw std::runtime_error(
      "ExplicitGpuFftStepper requires DIFFUSION3D_CUDA=ON");
#else
  settings.validate();
  species_ids_ = grid.species();
  if (species_ids_.empty())
    throw std::runtime_error("Grid must contain at least one species");

  settings_ = settings;
  cfg_nx_ = grid.nx_;
  cfg_ny_ = grid.ny_;
  cfg_nz_ = grid.nz_;
  cfg_h_ = grid.h_;

  const int fx = derive_fft_axis_len(cfg_nx_);
  const int fy = derive_fft_axis_len(cfg_ny_);
  const int fz = derive_fft_axis_len(cfg_nz_);

  for (auto id : species_ids_) {
    destroy_fft_scratch(id);
    fft_scratch_[id] =
        diffusion3d_fft_scratch_create(cfg_nx_, cfg_ny_, cfg_nz_, fx, fy, fz);

    const double D = settings.species_diffusivities.at(id);
    diffusion3d_fft_scratch_rebuild_operator(fft_scratch_[id], cfg_h_, D,
                                             tick_dt, settings.safety);
  }

  configured_ = true;
#endif
}

void ExplicitGpuFftStepper::advance_species_interval(
    MultiSpeciesFieldGrid &grid) {
#ifndef DIFFUSION3D_CUDA
  throw std::runtime_error(
      "ExplicitGpuFftStepper requires DIFFUSION3D_CUDA=ON");
#else
  if (!configured_)
    throw std::runtime_error(
        "Stepper not configured; call configure_species_interval() first");

  if (grid.nx_ != cfg_nx_ || grid.ny_ != cfg_ny_ || grid.nz_ != cfg_nz_ ||
      grid.h_ != cfg_h_)
    throw std::runtime_error(
        "Grid dimensions changed since configure_species_interval()");

  for (auto species_id : species_ids_) {
    ScalarFieldGrid &species_grid = grid.grid(species_id);

    diffusion3d_fft_scratch_apply(fft_scratch_.at(species_id),
                                  species_grid.data_);
  }
#endif
}

std::vector<SpeciesId> ExplicitGpuFftStepper::supported_species() const {
  return species_ids_;
}
