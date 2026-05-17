/**
 * @file multi_species_diffusion_engine.h
 * @brief Interface for multi-species diffusion engines in Diffusion3D.
 *
 * Defines the MultiSpeciesDiffusionEngine abstract base class and DiffusionAlgorithm enum.
 */

#ifndef DIFFUSION3D_MULTI_SPECIES_DIFFUSION_ENGINE_H
#define DIFFUSION3D_MULTI_SPECIES_DIFFUSION_ENGINE_H

// Multi-species diffusion engine interface
using SpeciesId = int; // Placeholder, replace with canonical SpeciesId later

enum class DiffusionAlgorithm
{
    ExplicitHeatEquation,
    // Add more algorithms as needed
};

class MultiSpeciesDiffusionEngine
{
public:
    virtual ~MultiSpeciesDiffusionEngine() = default;
    /**
     * @brief Configure the engine for a given species interval.
     */
    virtual void configure_species_interval(const MultiSpeciesFieldGrid &, const MultiSpeciesDiffusionSettings &) = 0;
    /**
     * @brief Advance all species by one interval.
     */
    virtual void advance_species_interval(MultiSpeciesFieldGrid &) = 0;
    /**
     * @brief Get the supported species IDs.
     */
    virtual std::vector<SpeciesId> supported_species() const = 0;
    /**
     * @brief Get the resolved algorithm type.
     */
    virtual DiffusionAlgorithm resolved_algorithm() const = 0;
};

#endif // DIFFUSION3D_MULTI_SPECIES_DIFFUSION_ENGINE_H
// Multi-species diffusion engine interface
#pragma once
#include <vector>
#include "multi_species_field_grid.h"
#include "multi_species_diffusion_settings.h"

using SpeciesId = int; // Placeholder, replace with canonical SpeciesId later

enum class DiffusionAlgorithm
{
    ExplicitHeatEquation
    // Add more algorithms as needed
};

class MultiSpeciesDiffusionEngine
{
public:
    virtual ~MultiSpeciesDiffusionEngine() = default;
    virtual void configure_species_interval(const MultiSpeciesFieldGrid &, const MultiSpeciesDiffusionSettings &) = 0;
    virtual void advance_species_interval(MultiSpeciesFieldGrid &) = 0;
    virtual std::vector<SpeciesId> supported_species() const = 0;
    virtual DiffusionAlgorithm resolved_algorithm() const = 0;
};
