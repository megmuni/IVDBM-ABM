#ifndef IVDBM_CHEMICAL_ENVIRONMENT_CONFIG_H
#define IVDBM_CHEMICAL_ENVIRONMENT_CONFIG_H

/**
 * @file chemical_environment_config.h
 * @brief JSON-backed settings for the IVDBM chemical environment.
 */

#include "species_id.h"

#include <map>
#include <string>
#include <vector>

/** One diffusing species entry from the config file. */
struct SpeciesConfigEntry
{
    SpeciesId id = -1;
    std::string name;
    double base_diffusivity_mm2_per_min = 0.0;
    int concentration_channel = -1;
    int diffused_channel = -1;
};

/**
 * @brief Full chemical environment definition loaded from JSON.
 *
 * Copy @c configFiles/chemical_environment.template.json to
 * @c configFiles/chemical_environment.json and edit; simulation code reads
 * only from the JSON path (see util::getChemicalEnvironmentConfigPath()).
 */
struct ChemicalEnvironmentConfig
{
    int schema_version = 0;
    std::string model;
    double tick_interval_minutes = 30.0;
    int channel_count = 0;
    int chemotaxis_channel = -1;
    std::string merge_chemotaxis_from_species;
    std::map<std::string, float> baseline_total_mass;
    std::vector<SpeciesConfigEntry> species;

    /** Total baseline mass for a species name; throws if missing. */
    float baseline_total_mass_for(const std::string &species_name) const;
};

/** Load and validate a chemical environment JSON file. */
ChemicalEnvironmentConfig load_chemical_environment_config(const std::string &path);

#endif
