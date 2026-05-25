#ifndef IVDBM_SPECIES_REGISTRY_H
#define IVDBM_SPECIES_REGISTRY_H

/**
 * @file species_registry.h
 * @brief Central registry of diffusing species and their physical parameters.
 *
 * Cytokine-specific names and chemAllocation channel indices live here (and in
 * factory helpers), not in PatchFieldDiffusion or diffusion3d_core.
 */

#include "species_id.h"

#include <map>
#include <string>
#include <vector>

struct SpeciesDescriptor
{
    SpeciesId id = -1;
    std::string name;
    /** Base diffusivity scale (mm²/min) before hydrogel swelling `Q` is applied. */
    double base_diffusivity = 0.0;
    /** chemAllocation index for present/concentration channel (e.g. pTNF). */
    int concentration_channel = -1;
    /** chemAllocation index for diffused output channel (e.g. dTNF). */
    int diffused_channel = -1;
};

class SpeciesRegistry
{
public:
    void set_swelling_ratio(double Q);
    void register_species(const SpeciesDescriptor &desc);
    void clear();

    const SpeciesDescriptor &descriptor(SpeciesId id) const;
    double diffusivity(SpeciesId id) const;
    std::vector<SpeciesId> diffusing_species() const;
    bool empty() const { return species_.empty(); }

    /** IVDBM baseline cytokines (TNF, TGF, IL-1β) with legacy diffusivity constants. */
    static SpeciesRegistry ivdbm_default(double swelling_ratio_Q);

private:
    double Q_ = 1.0;
    std::map<SpeciesId, SpeciesDescriptor> species_;
    std::vector<SpeciesId> order_;
};

#endif
