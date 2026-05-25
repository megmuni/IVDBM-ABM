#ifndef IVDBM_CHEMICAL_ENVIRONMENT_H
#define IVDBM_CHEMICAL_ENVIRONMENT_H

/**
 * @file chemical_environment.h
 * @brief Facade for multi-species chemistry: diffusion phase + tick merge
 * (Phase III).
 *
 * Owns no PDE logic (delegates to PatchFieldDiffusion). During III.2, buffers
 * remain legacy chemAllocation rows; WorldChem totals updated on merge.
 */

#include "../World/PatchFieldDiffusion.h"
#include "chemical_channel_views.h"
#include "chemotaxis_signal.h"
#include "species_registry.h"

#include <functional>
#include <map>

class Chemical;

class ChemicalEnvironment {
public:
  ChemicalEnvironment(int nx, int ny, int nz);

  /** Load IVDBM default cytokines and diffusivities for hydrogel swelling Q. */
  void load_ivdbm_default(double swelling_ratio_Q);

  void set_swelling_ratio(double Q);

  /** To link chemAllocation rows and chemotaxis channel to the environment. */
  void bind_legacy_grid(float **chem_allocation, int chemotaxis_channel_index);

  /**
   * Injects diffusion callback (typically
   * PatchFieldDiffusion::diffuse_all_species).
   */
  void set_diffusion_runner(std::function<void(double tick_dt_minutes)> runner);

  /** Stage 1: clear d*, then invoke diffusion runner. */
  void run_diffusion_phase(double tick_dt_minutes);

  /** Buffer map for the registered diffusion runner (p* / d* views). */
  std::map<SpeciesId, SpeciesDiffusionBuffers> diffusion_buffers() const;

  /**
   * Stage 4a: p* += d*, clear d*, refresh pcellgrad and WorldChem totals.
   * world_chem must remain linked to the same chemAllocation rows.
   */
  void merge_and_reset_secretion(Chemical &world_chem);

  SpeciesChannelViews channels(SpeciesId id) const;
  ChemotaxisSignal chemotaxis_signal() const;

  const SpeciesRegistry &registry() const { return registry_; }

  /** Agent API (III.6): local concentration p* at patch index. */
  float concentration_at(int patch_index, SpeciesId species) const;

  /** Agent API: accumulate into legacy d* (secretion + diffused pre-merge). */
  void accumulate_secretion(int patch_index, SpeciesId species, float delta);

  /** Agent API: read by legacy concentration channel index (e.g. pTNF). */
  float concentration_at_channel(int patch_index,
                                 int concentration_channel) const;

private:
  void zero_delta_channels();

  int nx_;
  int ny_;
  int nz_;
  int grid_size_;

  SpeciesRegistry registry_;
  std::function<void(double)> diffusion_runner_;

  float **chem_allocation_ = nullptr;
  int chemotaxis_channel_ = -1;
};

#endif
