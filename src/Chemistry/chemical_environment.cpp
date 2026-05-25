#include "chemical_environment.h"

#include "../FieldVariable/Usr_FieldVariables/Chemical.h"

#include <algorithm>
#include <stdexcept>

ChemicalEnvironment::ChemicalEnvironment(int nx, int ny, int nz)
    : nx_(nx), ny_(ny), nz_(nz), grid_size_(nx * ny * nz) {
  if (nx_ <= 0 || ny_ <= 0 || nz_ <= 0)
    throw std::invalid_argument(
        "ChemicalEnvironment: grid dimensions must be positive");
}

void ChemicalEnvironment::load_ivdbm_default(double swelling_ratio_Q) {
  registry_ = SpeciesRegistry::ivdbm_default(swelling_ratio_Q);
  if (registry_.empty())
    throw std::runtime_error(
        "ChemicalEnvironment: ivdbm_default produced empty registry");
}

void ChemicalEnvironment::set_swelling_ratio(double Q) {
  registry_.set_swelling_ratio(Q);
}

void ChemicalEnvironment::bind_legacy_grid(float **chem_allocation,
                                           int chemotaxis_channel_index) {
  if (chem_allocation == nullptr)
    throw std::invalid_argument("ChemicalEnvironment: chem_allocation is null");
  chem_allocation_ = chem_allocation;
  chemotaxis_channel_ = chemotaxis_channel_index;
}

void ChemicalEnvironment::set_diffusion_runner(
    std::function<void(double tick_dt_minutes)> runner) {
  diffusion_runner_ = std::move(runner);
}

void ChemicalEnvironment::zero_delta_channels() {
  for (SpeciesId id : registry_.diffusing_species()) {
    const int ch = registry_.descriptor(id).diffused_channel;
    std::fill(chem_allocation_[ch], chem_allocation_[ch] + grid_size_, 0.f);
  }
}

std::map<SpeciesId, SpeciesDiffusionBuffers>
ChemicalEnvironment::diffusion_buffers() const {
  std::map<SpeciesId, SpeciesDiffusionBuffers> buffers;
  for (SpeciesId id : registry_.diffusing_species()) {
    const SpeciesDescriptor &desc = registry_.descriptor(id);
    SpeciesDiffusionBuffers buf;
    buf.concentration = chem_allocation_[desc.concentration_channel];
    buf.diffused = chem_allocation_[desc.diffused_channel];
    buffers[id] = buf;
  }
  return buffers;
}

void ChemicalEnvironment::run_diffusion_phase(double tick_dt_minutes) {
  if (registry_.empty() || chem_allocation_ == nullptr)
    return;

  zero_delta_channels();
  if (diffusion_runner_)
    diffusion_runner_(tick_dt_minutes);
}

SpeciesChannelViews ChemicalEnvironment::channels(SpeciesId id) const {
  const SpeciesDescriptor &desc = registry_.descriptor(id);
  SpeciesChannelViews views;
  views.concentration = chem_allocation_[desc.concentration_channel];
  views.secretion_delta = chem_allocation_[desc.diffused_channel];
  views.diffused = views.secretion_delta;
  return views;
}

ChemotaxisSignal ChemicalEnvironment::chemotaxis_signal() const {
  ChemotaxisSignal signal;
  if (chem_allocation_ != nullptr && chemotaxis_channel_ >= 0)
    signal.data = chem_allocation_[chemotaxis_channel_];
  return signal;
}

float ChemicalEnvironment::concentration_at(int patch_index,
                                           SpeciesId species) const
{
    const SpeciesDescriptor &desc = registry_.descriptor(species);
    return chem_allocation_[desc.concentration_channel][patch_index];
}

void ChemicalEnvironment::accumulate_secretion(int patch_index,
                                               SpeciesId species, float delta)
{
    const SpeciesDescriptor &desc = registry_.descriptor(species);
    chem_allocation_[desc.diffused_channel][patch_index] += delta;
}

float ChemicalEnvironment::concentration_at_channel(
    int patch_index, int concentration_channel) const
{
    return chem_allocation_[concentration_channel][patch_index];
}

void ChemicalEnvironment::merge_and_reset_secretion(Chemical &world_chem) {
  if (registry_.empty() || chem_allocation_ == nullptr)
    return;

  world_chem.totalTNF = 0;
  world_chem.totalTGF = 0;
  world_chem.totalIL1beta = 0;

  float sum_tnf = 0, sum_tgf = 0, sum_il1 = 0;
  const std::vector<SpeciesId> diffusing = registry_.diffusing_species();

  for (int zi = 0; zi < nz_; ++zi) {
    for (int yi = 0; yi < ny_; ++yi) {
      for (int xi = 0; xi < nx_; ++xi) {
        const int in = xi + yi * nx_ + zi * nx_ * ny_;

        for (SpeciesId id : diffusing) {
          const SpeciesDescriptor &desc = registry_.descriptor(id);
          float *p = chem_allocation_[desc.concentration_channel];
          float *d = chem_allocation_[desc.diffused_channel];
#ifndef CALIBRATION
          p[in] = d[in] + p[in];
#else
          p[in] = d[in] + p[in] * 0.02f;
#endif
          d[in] = 0.f;
        }

        world_chem.pcellgrad[in] = world_chem.pTGF[in];

        sum_tnf += world_chem.pTNF[in];
        sum_tgf += world_chem.pTGF[in];
        sum_il1 += world_chem.pIL1beta[in];
      }
    }
  }

  world_chem.totalTNF += sum_tnf;
  world_chem.totalTGF += sum_tgf;
  world_chem.totalIL1beta += sum_il1;
}
