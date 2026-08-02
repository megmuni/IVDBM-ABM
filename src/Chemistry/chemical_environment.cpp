/**
 * @file chemical_environment.cpp
 * @brief Implements per-tick chemistry: storage, diffusion hook, merge, totals.
 */

#include "chemical_environment.h"

#include "../enums.h"

#include <algorithm>
#include <stdexcept>

ChemicalEnvironment::ChemicalEnvironment(int nx, int ny, int nz,
                                         double grid_spacing_mm)
    : nx_(nx), ny_(ny), nz_(nz), grid_size_(nx * ny * nz),
      grid_spacing_mm_(grid_spacing_mm) {
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

void ChemicalEnvironment::load_from_config(const std::string &config_path,
                                           double swelling_ratio_Q,
                                           double stiffness_E) {
  config_ = load_chemical_environment_config(config_path);
  registry_ =
      SpeciesRegistry::from_config(config_, swelling_ratio_Q, stiffness_E);
  merge_chemotaxis_species_ =
      species_id_by_name(config_.merge_chemotaxis_from_species);
  sync_diffusion_registry();
}

void ChemicalEnvironment::set_swelling_ratio(double Q) {
  registry_.set_swelling_ratio(Q);
}

void ChemicalEnvironment::set_stiffness(double E) {
  registry_.set_stiffness(E);
}

float ChemicalEnvironment::baseline_total_mass_for(
    const std::string &species_name) const {
  return config_.baseline_total_mass_for(species_name);
}

const char *ChemicalEnvironment::diffusion_algorithm_label() const {
  if (diffusion_runner_)
    return "custom diffusion runner";
  if (!patch_diffusion_)
    return "diffusion disabled (grid spacing not set)";
  return patch_diffusion_->effective_algorithm_label();
}

DiffusionAlgorithm ChemicalEnvironment::diffusion_algorithm() const {
  if (patch_diffusion_)
    return patch_diffusion_->configured_algorithm();
  return patch_field_default_diffusion_algorithm();
}

void ChemicalEnvironment::set_diffusion_algorithm(DiffusionAlgorithm algo) {
  if (patch_diffusion_)
    patch_diffusion_->set_diffusion_algorithm(algo);
}

void ChemicalEnvironment::allocate_channels_from_config() {
  if (config_.channel_count <= 0)
    throw std::runtime_error("ChemicalEnvironment: load_from_config before "
                             "allocate_channels_from_config");
  allocate_channel_storage(config_.channel_count, config_.chemotaxis_channel);
}

SpeciesId
ChemicalEnvironment::species_id_by_name(const std::string &name) const {
  for (const SpeciesConfigEntry &s : config_.species) {
    if (s.name == name)
      return s.id;
  }
  throw std::out_of_range("ChemicalEnvironment: unknown species name " + name);
}

int ChemicalEnvironment::concentration_channel_for(
    const std::string &name) const {
  return registry_.descriptor(species_id_by_name(name)).concentration_channel;
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
  total_o2_ = 0.f;

  if (channel_data_.empty() || registry_.empty())
    return;

  const int tnf_ch = concentration_channel_for("TNF");
  const int tgf_ch = concentration_channel_for("TGF");
  const int il1_ch = concentration_channel_for("IL1beta");
  const int o2_ch = concentration_channel_for("o2");

  for (int i = 0; i < grid_size_; ++i) {
    total_tnf_ += channel_row(tnf_ch)[i];
    total_tgf_ += channel_row(tgf_ch)[i];
    total_il1beta_ += channel_row(il1_ch)[i];
    total_o2_ += channel_row(o2_ch)[i];
  }
}

std::map<SpeciesId, SpeciesDiffusionBuffers>
ChemicalEnvironment::diffusion_buffers() const {
  std::map<SpeciesId, SpeciesDiffusionBuffers> buffers;
  for (SpeciesId id : registry_.diffusing_species()) {
    const SpeciesDescriptor &desc = registry_.descriptor(id);
    SpeciesDiffusionBuffers buf;
    buf.concentration =
        const_cast<float *>(channel_row(desc.concentration_channel));
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
  views.secretion_delta =
      const_cast<float *>(channel_row(desc.diffused_channel));
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
  total_o2_ = 0.f;

  const std::vector<SpeciesId> diffusing = registry_.diffusing_species();
  const int tnf_ch = concentration_channel_for("TNF");
  const int tgf_ch = concentration_channel_for("TGF");
  const int il1_ch = concentration_channel_for("IL1beta");
  const int o2_ch = concentration_channel_for("o2");
  const int chemo_src =
      registry_.descriptor(merge_chemotaxis_species_).concentration_channel;

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
          p[in] = std::max(p[in], 0.f);
          d[in] = 0.f;
        }

        if (chemotaxis_channel_ >= 0) {
          float chemo = channel_row(chemo_src)[in];
          channel_row(chemotaxis_channel_)[in] = std::max(chemo, 0.f);
        }

        total_tnf_ += channel_row(tnf_ch)[in];
        total_tgf_ += channel_row(tgf_ch)[in];
        total_il1beta_ += channel_row(il1_ch)[in];
        total_o2_ += channel_row(o2_ch)[in];
      }
    }
  }
}
