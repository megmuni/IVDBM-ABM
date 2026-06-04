#ifndef IVDBM_CHEMICAL_ENVIRONMENT_H
#define IVDBM_CHEMICAL_ENVIRONMENT_H

/**
 * @file chemical_environment.h
 * @brief Coordinates cytokine storage and per-tick updates on the patch grid.
 *
 * Biological parameters (species, diffusivity, baselines, tick length) are
 * loaded from JSON - see configFiles/chemical_environment.template.json.
 */

#include "../Diffusion3D/core/diffusion_algorithm.h"
#include "chemical_channel_views.h"
#include "chemical_environment_config.h"
#include "chemotaxis_signal.h"
#include "diffusion_algorithm_name.h"
#include "patch_field_diffusion.h"
#include "species_registry.h"

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

class Chemical;

/**
 * @brief Per-tick chemistry for BMWorld: storage, diffusion, merge, totals.
 */
class ChemicalEnvironment {
public:
  /**
   * @param grid_spacing_mm Patch length in mm; pass a value > 0 to enable
   * diffusion. Tests may omit this (default 0) when only merge/storage is
   * needed.
   */
  ChemicalEnvironment(int nx, int ny, int nz, double grid_spacing_mm = 0.0);

  /**
   * @brief Load species, channels, and baselines from JSON.
   * @param config_path Path to chemical_environment.json (copy from template).
   * @param swelling_ratio_Q Hydrogel swelling ratio applied to diffusivity.
   */
  void load_from_config(const std::string &config_path,
                        double swelling_ratio_Q);

  /** Update effective diffusivity when hydrogel swelling ratio @p Q changes. */
  void set_swelling_ratio(double Q);

  /** Settings last passed to load_from_config. */
  const ChemicalEnvironmentConfig &configuration() const { return config_; }

  /** Macro tick length in minutes (from config). */
  double tick_interval_minutes() const { return config_.tick_interval_minutes; }

  /**
   * @brief Human-readable diffusion algorithm (or why diffusion is off).
   * Shows the effective solver (resolved backend), not only the requested one.
   */
  const char *diffusion_algorithm_label() const;

  /** Requested Diffusion3D algorithm for patch-field diffusion. */
  DiffusionAlgorithm diffusion_algorithm() const;

  /**
   * @brief Select the patch-field diffusion algorithm (ignored if a custom
   *        diffusion runner is set via set_diffusion_runner).
   */
  void set_diffusion_algorithm(DiffusionAlgorithm algo);

  /** Total baseline mass for a species (from config), e.g. @c "TNF". */
  float baseline_total_mass_for(const std::string &species_name) const;

  /**
   * @brief Allocate channel grids (usually call after load_from_config).
   * Uses channel count and chemotaxis index from the config file.
   */
  void allocate_channels_from_config();

  /** @deprecated Prefer allocate_channels_from_config after load_from_config.
   */
  void allocate_channel_storage(int channel_count,
                                int chemotaxis_channel_index);

  int channel_count() const { return static_cast<int>(channel_data_.size()); }

  void set_concentration(int patch_index, SpeciesId species, float value);

  void clear_delta_channels();

  void update_chemotaxis_from_species(SpeciesId species);

  void recompute_world_totals();

  void set_diffusion_runner(std::function<void(double tick_dt_minutes)> runner);

  void run_diffusion_phase(double tick_dt_minutes);

  std::map<SpeciesId, SpeciesDiffusionBuffers> diffusion_buffers() const;

  void merge_and_reset_secretion();

  void copy_totals_to(Chemical &world_chem) const;

  SpeciesChannelViews channels(SpeciesId id) const;
  ChemotaxisSignal chemotaxis_signal() const;

  const SpeciesRegistry &registry() const { return registry_; }

  float concentration_at(int patch_index, SpeciesId species) const;

  void accumulate_secretion(int patch_index, SpeciesId species, float delta);

  float concentration_at_channel(int patch_index,
                                 int concentration_channel) const;

  float chemotaxis_at(int patch_index) const;

  const float *channel_grid(int channel_index) const;

  int grid_nx() const { return nx_; }
  int grid_ny() const { return ny_; }
  int grid_nz() const { return nz_; }
  int grid_size() const { return grid_size_; }

  /** Patch spacing in mm (from construction); 0 if diffusion was not enabled. */
  double grid_spacing_mm() const { return grid_spacing_mm_; }

  /** Chemotaxis channel index from config, or -1 before allocation. */
  int chemotaxis_channel_index() const { return chemotaxis_channel_; }

  float total_tnf() const { return total_tnf_; }
  float total_tgf() const { return total_tgf_; }
  float total_il1beta() const { return total_il1beta_; }

private:
  void sync_diffusion_registry();
  SpeciesId species_id_by_name(const std::string &name) const;
  int concentration_channel_for(const std::string &name) const;

  float *channel_row(int channel);
  const float *channel_row(int channel) const;

  int nx_;
  int ny_;
  int nz_;
  int grid_size_;
  double grid_spacing_mm_ = 0.0;

  ChemicalEnvironmentConfig config_;
  SpeciesRegistry registry_;
  SpeciesId merge_chemotaxis_species_ = -1;
  std::unique_ptr<PatchFieldDiffusion> patch_diffusion_;
  std::function<void(double)> diffusion_runner_;

  std::vector<std::vector<float>> channel_data_;
  int chemotaxis_channel_ = -1;

  float total_tnf_ = 0.f;
  float total_tgf_ = 0.f;
  float total_il1beta_ = 0.f;
};

#endif
