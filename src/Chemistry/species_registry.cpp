#include "species_registry.h"

#include "../enums.h"

#include <stdexcept>

void SpeciesRegistry::set_swelling_ratio(double Q)
{
    if (Q <= 0.0)
        throw std::invalid_argument("SpeciesRegistry: swelling ratio Q must be > 0");
    Q_ = Q;
}

void SpeciesRegistry::set_stiffness(double E)
{
    if (E <= 0.0)
        throw std::invalid_argument("SpeciesRegistry: stiffness E must be > 0");
    E_ = E;
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
    //TO-DO: modify to check for which type of effective D calculation (i.e., using Q or using E)
    //is desired for each chemical. For the 3 current cytokines, Q should be used. For O2, E should
    //be used.
}

std::vector<SpeciesId> SpeciesRegistry::diffusing_species() const
{
    return order_;
}

SpeciesRegistry SpeciesRegistry::from_config(const ChemicalEnvironmentConfig &cfg,
                                             double swelling_ratio_Q, double stiffness_E)
{
    SpeciesRegistry registry;
    registry.set_swelling_ratio(swelling_ratio_Q);
    registry.set_stiffness(stiffness_E);

    for (const SpeciesConfigEntry &s : cfg.species)
    {
        SpeciesDescriptor desc;
        desc.id = s.id;
        desc.name = s.name;
        desc.base_diffusivity = s.base_diffusivity_mm2_per_min;
        desc.concentration_channel = s.concentration_channel;
        desc.diffused_channel = s.diffused_channel;
        registry.register_species(desc);
    }

    return registry;
}
