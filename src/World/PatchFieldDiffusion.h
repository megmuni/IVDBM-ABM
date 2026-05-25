#ifndef IVDBM_PATCH_FIELD_DIFFUSION_H
#define IVDBM_PATCH_FIELD_DIFFUSION_H

/**
 * @file PatchFieldDiffusion.h
 * @brief World adapter: ABM float buffers ↔ Diffusion3D CPU engine (Phase II).
 *
 * Interface only in II.0–II.1; implementation lands in II.3. No CUDA/GPU paths.
 */

#include "../Chemistry/species_id.h"
#include "../Chemistry/species_registry.h"

#include <cstddef>
#include <map>
#include <memory>

/** Per-species ABM buffer handles (size = nx * ny * nz). */
struct SpeciesDiffusionBuffers
{
    const float *concentration = nullptr;
    float *diffused = nullptr;
};

class PatchFieldDiffusion
{
public:
    PatchFieldDiffusion(int nx, int ny, int nz, double grid_spacing);
    ~PatchFieldDiffusion();

    PatchFieldDiffusion(const PatchFieldDiffusion &) = delete;
    PatchFieldDiffusion &operator=(const PatchFieldDiffusion &) = delete;

    void set_registry(const SpeciesRegistry &registry);

    /**
     * Advance all diffusing species over one macro tick (`tick_dt`, same units as BMWorld).
     * Reads `concentration`, writes incremental diffused field into `diffused`.
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
