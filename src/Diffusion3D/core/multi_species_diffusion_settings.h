/**
 * @file multi_species_diffusion_settings.h
 * @brief Multi-species configuration with per-species diffusivities for Diffusion3D.
 *
 * Defines MultiSpeciesDiffusionSettings for species-specific diffusion parameters.
 */

#ifndef DIFFUSION3D_MULTI_SPECIES_DIFFUSION_SETTINGS_H
#define DIFFUSION3D_MULTI_SPECIES_DIFFUSION_SETTINGS_H

// Multi-species configuration with per-species diffusivities
#include <map>
#include <string>
#include <vector>
#include <stdexcept>

using SpeciesId = int; // Placeholder, replace with canonical SpeciesId later

class MultiSpeciesDiffusionSettings
{
public:
    std::map<SpeciesId, double> species_diffusivities_;
    std::vector<std::string> species_names_;

    void validate() const
    {
        if (species_diffusivities_.empty())
            throw std::runtime_error("No species diffusivities provided");
        for (const auto &[id, d] : species_diffusivities_)
        {
            if (d <= 0.0)
                throw std::runtime_error("Diffusivity must be positive");
        }
    }
};

#endif // DIFFUSION3D_MULTI_SPECIES_DIFFUSION_SETTINGS_H
