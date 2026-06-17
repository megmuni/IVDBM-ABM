#ifndef IVDBM_PATCH_FIELD_DIFFUSION_H
#define IVDBM_PATCH_FIELD_DIFFUSION_H

/**
 * @file patch_field_diffusion.h
 * @brief Connects ABM patch buffers to the Diffusion3D CPU solver.
 */

#include "../Diffusion3D/core/diffusion_algorithm.h"
#include "diffusion_algorithm_name.h"
#include "species_id.h"
#include "species_registry.h"

#include <cstddef>
#include <map>
#include <memory>

/** Per-species buffer handles on the patch grid (length nx × ny × nz). */
struct SpeciesDiffusionBuffers
{
    const float *concentration = nullptr;
    float *diffused = nullptr;
};

/**
 * @brief Runs one macro-tick of 3D diffusion for all registered species.
 */
class PatchFieldDiffusion
{
public:
    PatchFieldDiffusion(int nx, int ny, int nz, double grid_spacing);
    ~PatchFieldDiffusion();

    PatchFieldDiffusion(const PatchFieldDiffusion &) = delete;
    PatchFieldDiffusion &operator=(const PatchFieldDiffusion &) = delete;

    void set_registry(const SpeciesRegistry &registry);

    /**
     * @brief Advance all diffusing species over one tick.
     * @param buffers Concentration (read) and diffused row (write increment).
     * @param tick_dt Macro tick length in minutes (same as BMWorld).
     */
    void diffuse_all_species(const std::map<SpeciesId, SpeciesDiffusionBuffers> &buffers,
                             double tick_dt);

    /** Requested algorithm (see patch_field_default_diffusion_algorithm() for initial value). */
    void set_diffusion_algorithm(DiffusionAlgorithm algo);
    DiffusionAlgorithm configured_algorithm() const;

    /** Label for configured_algorithm() (see diffusion_algorithm_name.h). */
    const char *configured_algorithm_label() const;

    /** Label for the solver actually used (may differ if CUDA is off). */
    const char *effective_algorithm_label() const;

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    double h() const { return h_; }

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;

    int nx_;
    int ny_;
    int nz_;
    double h_;
    DiffusionAlgorithm algorithm_ = patch_field_default_diffusion_algorithm();
    const SpeciesRegistry *registry_ = nullptr;
};

#endif
