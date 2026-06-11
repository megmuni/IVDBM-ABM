#ifndef IVDBM_SPECIES_REGISTRY_H
#define IVDBM_SPECIES_REGISTRY_H

/**
 * @file species_registry.h
 * @brief Lists diffusing species and maps them to grid channels and diffusivity.
 *
 * Keeps cytokine names and channel indices out of PatchFieldDiffusion and
 * Diffusion3D core code.
 */

#include "chemical_environment_config.h"
#include "species_id.h"

#include <map>
#include <string>
#include <vector>

/** Metadata for one diffusing species on the patch grid. */
struct SpeciesDescriptor
{
    SpeciesId id = -1;
    std::string name;
    /** Base diffusivity (mm²/min) before swelling ratio Q is applied. */
    double base_diffusivity = 0.0;
    /** Grid channel for stored concentration (e.g. pTNF). */
    int concentration_channel = -1;
    /** Grid channel for per-tick delta / diffusion output (e.g. dTNF). */
    int diffused_channel = -1;
};

/**
 * @brief Registry of species, effective D, and channel indices.
 */
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

    /** Build registry from a validated @ref ChemicalEnvironmentConfig. */
    static SpeciesRegistry from_config(const ChemicalEnvironmentConfig &cfg,
                                       double swelling_ratio_Q);

private:
    double Q_ = 1.0;
    std::map<SpeciesId, SpeciesDescriptor> species_;
    std::vector<SpeciesId> order_;
};

#endif
