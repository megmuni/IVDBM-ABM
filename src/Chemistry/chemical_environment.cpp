/**
 * @file chemical_environment.cpp
 * @brief Implements per-tick chemistry: storage, diffusion hook, merge, totals.
 */

#include "chemical_environment.h"

#include "../enums.h"
#include "../FieldVariable/Usr_FieldVariables/Chemical.h"

#include <algorithm>
#include <stdexcept>

ChemicalEnvironment::ChemicalEnvironment(int nx, int ny, int nz,
                                           double grid_spacing_mm)
    : nx_(nx), ny_(ny), nz_(nz), grid_size_(nx * ny * nz) {
  if (nx_ <= 0 || ny_ <= 0 || nz_ <= 0)
    throw std::invalid_argument(
        "ChemicalEnvironment: grid dimensions must be positive");
  if (grid_spacing_mm > 0.0)
    patch_diffusion_ =
        std::make_unique<PatchFieldDiffusion>(nx_, ny_, nz_, grid_spacing_mm);
}

void ChemicalEnvironment::sync_diffusion_registry() {
  if (patch_diffusion_ && !registry_.empty())
    patch_diffusion_->set_registry(registry_);
}

void ChemicalEnvironment::load_ivdbm_default(double swelling_ratio_Q) {
  registry_ = SpeciesRegistry::ivdbm_default(swelling_ratio_Q);
  if (registry_.empty())
    throw std::runtime_error(
        "ChemicalEnvironment: ivdbm_default produced empty registry");
  sync_diffusion_registry();
}

void ChemicalEnvironment::set_swelling_ratio(double Q) {
  registry_.set_swelling_ratio(Q);
}

void ChemicalEnvironment::allocate_channel_storage(int channel_count,
                                                    int chemotaxis_channel) {
  if (channel_count <= 0)
    throw std::invalid_argument(
        "ChemicalEnvironment: channel_count must be positive");
  if (chemotaxis_channel < 0 || chemotaxis_channel >= channel_count)
    throw std::invalid_argument(
        "ChemicalEnvironment: invalid chemotaxis channel index");

  channel_data_.assign(static_cast<size_t>(channel_count),
                       std::vector<float>(grid_size_, 0.f));
  chemotaxis_channel_ = chemotaxis_channel;
}

float *ChemicalEnvironment::channel_row(int channel) {
  return channel_data_.at(static_cast<size_t>(channel)).data();
}

const float *ChemicalEnvironment::channel_row(int channel) const {
  return channel_data_.at(static_cast<size_t>(channel)).data();
}

void ChemicalEnvironment::set_diffusion_runner(
    std::function<void(double tick_dt_minutes)> runner) {
  diffusion_runner_ = std::move(runner);
}

void ChemicalEnvironment::set_concentration(int patch_index, SpeciesId species,
                                            float value) {
  const SpeciesDescriptor &desc = registry_.descriptor(species);
  channel_row(desc.concentration_channel)[patch_index] = value;
}

void ChemicalEnvironment::clear_delta_channels() {
  for (SpeciesId id : registry_.diffusing_species()) {
    const int ch = registry_.descriptor(id).diffused_channel;
    std::fill(channel_row(ch), channel_row(ch) + grid_size_, 0.f);
  }
}

void ChemicalEnvironment::update_chemotaxis_from_species(SpeciesId species) {
  if (channel_data_.empty() || chemotaxis_channel_ < 0)
    return;

  const int conc_ch = registry_.descriptor(species).concentration_channel;
  float *chemo = channel_row(chemotaxis_channel_);
  const float *conc = channel_row(conc_ch);
  for (int i = 0; i < grid_size_; ++i)
    chemo[i] = conc[i];
}

void ChemicalEnvironment::recompute_world_totals() {
  total_tnf_ = 0.f;
  total_tgf_ = 0.f;
  total_il1beta_ = 0.f;

  if (channel_data_.empty() || registry_.empty())
    return;

  const int tnf_ch = registry_.descriptor(TNF).concentration_channel;
  const int tgf_ch = registry_.descriptor(TGF).concentration_channel;
  const int il1_ch = registry_.descriptor(IL1beta).concentration_channel;

  for (int i = 0; i < grid_size_; ++i) {
    total_tnf_ += channel_row(tnf_ch)[i];
    total_tgf_ += channel_row(tgf_ch)[i];
    total_il1beta_ += channel_row(il1_ch)[i];
  }
}

std::map<SpeciesId, SpeciesDiffusionBuffers>
ChemicalEnvironment::diffusion_buffers() const {
  std::map<SpeciesId, SpeciesDiffusionBuffers> buffers;
  for (SpeciesId id : registry_.diffusing_species()) {
    const SpeciesDescriptor &desc = registry_.descriptor(id);
    SpeciesDiffusionBuffers buf;
    buf.concentration = const_cast<float *>(channel_row(desc.concentration_channel));
    buf.diffused = const_cast<float *>(channel_row(desc.diffused_channel));
    buffers[id] = buf;
  }
  return buffers;
}

void ChemicalEnvironment::run_diffusion_phase(double tick_dt_minutes) {
  if (registry_.empty() || channel_data_.empty())
    return;

  clear_delta_channels();
  if (diffusion_runner_) {
    diffusion_runner_(tick_dt_minutes);
    return;
  }
  if (patch_diffusion_)
    patch_diffusion_->diffuse_all_species(diffusion_buffers(), tick_dt_minutes);
}

SpeciesChannelViews ChemicalEnvironment::channels(SpeciesId id) const {
  const SpeciesDescriptor &desc = registry_.descriptor(id);
  SpeciesChannelViews views;
  views.concentration =
      const_cast<float *>(channel_row(desc.concentration_channel));
  views.secretion_delta = const_cast<float *>(channel_row(desc.diffused_channel));
  views.diffused = views.secretion_delta;
  return views;
}

ChemotaxisSignal ChemicalEnvironment::chemotaxis_signal() const {
  ChemotaxisSignal signal;
  if (!channel_data_.empty() && chemotaxis_channel_ >= 0)
    signal.data = const_cast<float *>(channel_row(chemotaxis_channel_));
  return signal;
}

float ChemicalEnvironment::concentration_at(int patch_index,
                                           SpeciesId species) const {
  const SpeciesDescriptor &desc = registry_.descriptor(species);
  return channel_row(desc.concentration_channel)[patch_index];
}

void ChemicalEnvironment::accumulate_secretion(int patch_index,
                                               SpeciesId species, float delta) {
  const SpeciesDescriptor &desc = registry_.descriptor(species);
  channel_row(desc.diffused_channel)[patch_index] += delta;
}

float ChemicalEnvironment::concentration_at_channel(
    int patch_index, int concentration_channel) const {
  return channel_row(concentration_channel)[patch_index];
}

float ChemicalEnvironment::chemotaxis_at(int patch_index) const {
  if (channel_data_.empty() || chemotaxis_channel_ < 0)
    return 0.f;
  return channel_row(chemotaxis_channel_)[patch_index];
}

const float *ChemicalEnvironment::channel_grid(int channel_index) const {
  return channel_row(channel_index);
}

void ChemicalEnvironment::merge_and_reset_secretion() {
  if (registry_.empty() || channel_data_.empty())
    return;

  total_tnf_ = 0.f;
  total_tgf_ = 0.f;
  total_il1beta_ = 0.f;

  const std::vector<SpeciesId> diffusing = registry_.diffusing_species();
  const int tnf_ch = registry_.descriptor(TNF).concentration_channel;
  const int tgf_ch = registry_.descriptor(TGF).concentration_channel;
  const int il1_ch = registry_.descriptor(IL1beta).concentration_channel;

  for (int zi = 0; zi < nz_; ++zi) {
    for (int yi = 0; yi < ny_; ++yi) {
      for (int xi = 0; xi < nx_; ++xi) {
        const int in = xi + yi * nx_ + zi * nx_ * ny_;

        for (SpeciesId id : diffusing) {
          const SpeciesDescriptor &desc = registry_.descriptor(id);
          float *p = channel_row(desc.concentration_channel);
          float *d = channel_row(desc.diffused_channel);
#ifndef CALIBRATION
          p[in] = d[in] + p[in];
#else
          p[in] = d[in] + p[in] * 0.02f;
#endif
          d[in] = 0.f;
        }

        if (chemotaxis_channel_ >= 0)
          channel_row(chemotaxis_channel_)[in] = channel_row(tgf_ch)[in];

        total_tnf_ += channel_row(tnf_ch)[in];
        total_tgf_ += channel_row(tgf_ch)[in];
        total_il1beta_ += channel_row(il1_ch)[in];
      }
    }
  }
}

void ChemicalEnvironment::copy_totals_to(Chemical &world_chem) const {
  world_chem.totalTNF = total_tnf_;
  world_chem.totalTGF = total_tgf_;
  world_chem.totalIL1beta = total_il1beta_;
}
