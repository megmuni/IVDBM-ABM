#ifndef IVDBM_CHEMICAL_ENVIRONMENT_H
#define IVDBM_CHEMICAL_ENVIRONMENT_H

/**
 * @file chemical_environment.h
 * @brief Coordinates cytokine storage and per-tick updates on the patch grid.
 *
 * Owns channel grids, the species registry, and (when grid spacing is given)
 * PatchFieldDiffusion for the default diffusion step each tick.
 */

#include "chemical_channel_views.h"
#include "chemotaxis_signal.h"
#include "patch_field_diffusion.h"
#include "species_registry.h"

#include <functional>
#include <map>
#include <memory>
#include <vector>

class Chemical;

/**
 * @brief Per-tick chemistry for BMWorld: storage, diffusion, merge, totals.
 * Is a singleton for the entire world.
 */
class ChemicalEnvironment {
public:
  /**
   * @param grid_spacing_mm Patch length in mm; pass a value > 0 to enable
   * diffusion. Tests may omit this (default 0) when only merge/storage is
   * needed.
   */
  ChemicalEnvironment(int nx, int ny, int nz, double grid_spacing_mm = 0.0);

  /** Register TNF, TGF, IL-1β with IVDBM default diffusivities scaled by @p
   * swelling_ratio_Q. */
  void load_ivdbm_default(double swelling_ratio_Q);

  /** Update effective diffusivity when hydrogel swelling ratio @p Q changes. */
  void set_swelling_ratio(double Q);

  /**
   * @brief Allocate one float grid per channel (size nx×ny×nz).
   * @param channel_count Total channels (concentration, delta, chemotaxis, …).
   * @param chemotaxis_channel_index Channel index cells use for chemotaxis
   * (e.g. pcellgrad).
   */
  void allocate_channel_storage(int channel_count,
                                int chemotaxis_channel_index);

  /** Number of allocated channels (zero until allocate_channel_storage). */
  int channel_count() const { return static_cast<int>(channel_data_.size()); }

  /** Set stored concentration for one species at a patch index. */
  void set_concentration(int patch_index, SpeciesId species, float value);

  /** Zero all per-tick delta (d*) channels before diffusion. */
  void clear_delta_channels();

  /** Copy one species’ concentration into the chemotaxis channel (e.g. TGF for
   * cell guidance). */
  void update_chemotaxis_from_species(SpeciesId species);

  /** Sum concentrations over the grid into internal totalTNF/TGF/IL1β. */
  void recompute_world_totals();

  /**
   * @brief Optional override for the diffusion step (tests or experiments).
   * When not set, run_diffusion_phase uses
   * PatchFieldDiffusion::diffuse_all_species.
   */
  void set_diffusion_runner(std::function<void(double tick_dt_minutes)> runner);

  /** Clear deltas, then run diffusion (built-in or custom runner). */
  void run_diffusion_phase(double tick_dt_minutes);

  /** Non-owning pointers into concentration and delta rows for the diffusion
   * adapter. */
  std::map<SpeciesId, SpeciesDiffusionBuffers> diffusion_buffers() const;

  /** Commit p* += d*, clear d*, refresh chemotaxis from TGF, update totals. */
  void merge_and_reset_secretion();

  /** Write integrated masses into @p world_chem (totals only; for CSV and
   * calibration). */
  void copy_totals_to(Chemical &world_chem) const;

  /** Non-owning views of concentration and delta rows for one species. */
  SpeciesChannelViews channels(SpeciesId id) const;

  /** Non-owning view of the chemotaxis channel. */
  ChemotaxisSignal chemotaxis_signal() const;

  const SpeciesRegistry &registry() const { return registry_; }

  /** Concentration of @p species at @p patch_index. */
  float concentration_at(int patch_index, SpeciesId species) const;

  /** Add @p delta to the per-tick secretion buffer for @p species at @p
   * patch_index. */
  void accumulate_secretion(int patch_index, SpeciesId species, float delta);

  /** Concentration using a raw channel index (see enums.h). */
  float concentration_at_channel(int patch_index,
                                 int concentration_channel) const;

  /** Chemotaxis field value at @p patch_index. */
  float chemotaxis_at(int patch_index) const;

  /** Row-major values for one channel (VTK export, debugging). */
  const float *channel_grid(int channel_index) const;

  float total_tnf() const { return total_tnf_; }
  float total_tgf() const { return total_tgf_; }
  float total_il1beta() const { return total_il1beta_; }

private:
  void sync_diffusion_registry();

  float *channel_row(int channel);
  const float *channel_row(int channel) const;

  int nx_;
  int ny_;
  int nz_;
  int grid_size_;

  SpeciesRegistry registry_;
  std::unique_ptr<PatchFieldDiffusion> patch_diffusion_;
  std::function<void(double)> diffusion_runner_;

  std::vector<std::vector<float>> channel_data_;
  int chemotaxis_channel_ = -1;

  float total_tnf_ = 0.f;
  float total_tgf_ = 0.f;
  float total_il1beta_ = 0.f;
};

#endif
