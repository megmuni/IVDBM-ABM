#include "species_registry.h"

#include "../enums.h"

#include <stdexcept>

void SpeciesRegistry::set_swelling_ratio(double Q)
{
    if (Q <= 0.0)
        throw std::invalid_argument("SpeciesRegistry: swelling ratio Q must be > 0");
    Q_ = Q;
}

void SpeciesRegistry::register_species(const SpeciesDescriptor &desc)
{
    if (desc.id < 0)
        throw std::invalid_argument("SpeciesRegistry: species id must be >= 0");
    if (desc.base_diffusivity <= 0.0)
        throw std::invalid_argument("SpeciesRegistry: base_diffusivity must be > 0");
    if (desc.concentration_channel < 0 || desc.diffused_channel < 0)
        throw std::invalid_argument("SpeciesRegistry: channel indices must be >= 0");

    if (species_.find(desc.id) == species_.end())
        order_.push_back(desc.id);
    species_[desc.id] = desc;
}

void SpeciesRegistry::clear()
{
    species_.clear();
    order_.clear();
}

const SpeciesDescriptor &SpeciesRegistry::descriptor(SpeciesId id) const
{
    const auto it = species_.find(id);
    if (it == species_.end())
        throw std::out_of_range("SpeciesRegistry: unknown species id");
    return it->second;
}

double SpeciesRegistry::diffusivity(SpeciesId id) const
{
    return descriptor(id).base_diffusivity * Q_;
}

std::vector<SpeciesId> SpeciesRegistry::diffusing_species() const
{
    return order_;
}

SpeciesRegistry SpeciesRegistry::ivdbm_default(double swelling_ratio_Q)
{
    SpeciesRegistry registry;
    registry.set_swelling_ratio(swelling_ratio_Q);

    // Matches legacy diffuseCytokines() coefficients: D = base * Q (mm²/min).
    const struct
    {
        SpeciesId id;
        const char *name;
        double base_D;
        int p_channel;
        int d_channel;
    } kIvdbm[] = {
        {TNF, "TNF", 0.0018 * 0.1, pTNF, dTNF},
        {TGF, "TGF", 0.00156 * 0.1, pTGF, dTGF},
        {IL1beta, "IL1beta", 0.0018 * 0.1, pIL1beta, dIL1beta},
    };

    for (const auto &s : kIvdbm)
    {
        SpeciesDescriptor desc;
        desc.id = s.id;
        desc.name = s.name;
        desc.base_diffusivity = s.base_D;
        desc.concentration_channel = s.p_channel;
        desc.diffused_channel = s.d_channel;
        registry.register_species(desc);
    }

    return registry;
}
