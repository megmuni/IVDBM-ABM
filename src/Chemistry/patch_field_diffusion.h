#ifndef IVDBM_PATCH_FIELD_DIFFUSION_H
#define IVDBM_PATCH_FIELD_DIFFUSION_H

/**
 * @file patch_field_diffusion.h
 * @brief Connects ABM patch buffers to the Diffusion3D CPU solver.
 */

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
    const SpeciesRegistry *registry_ = nullptr;
};

#endif
