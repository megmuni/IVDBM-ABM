#include "chemical_environment_config.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <stdexcept>

namespace
{

using json = nlohmann::json;

DiffusivityModelConfig parse_diffusivity_model(const json &entry,
                                               const std::string &species_name)
{
    DiffusivityModelConfig model;
    if (!entry.contains("diffusivity_model"))
        return model;

    const json &dm = entry.at("diffusivity_model");
    if (dm.is_string())
    {
        const std::string type = dm.get<std::string>();
        if (type != "swelling_ratio")
            throw std::invalid_argument(
                "chemical environment config: unknown diffusivity_model for " +
                species_name);
        return model;
    }

    if (!dm.is_object())
        throw std::invalid_argument(
            "chemical environment config: diffusivity_model must be string or "
            "object for " +
            species_name);

    const std::string type = dm.at("type").get<std::string>();
    if (type == "swelling_ratio")
        return model;

    if (type == "logarithmic_stiffness")
    {
        model.type = DiffusivityModelType::LogarithmicStiffness;
        model.slope = dm.at("slope").get<double>();
        model.intercept = dm.at("intercept").get<double>();
        return model;
    }

    throw std::invalid_argument(
        "chemical environment config: unknown diffusivity_model type for " +
        species_name);
}

} // namespace

float ChemicalEnvironmentConfig::baseline_total_mass_for(
    const std::string &species_name) const
{
    const auto it = baseline_total_mass.find(species_name);
    if (it == baseline_total_mass.end())
        throw std::out_of_range(
            "ChemicalEnvironmentConfig: unknown baseline species " +
            species_name);
    return it->second;
}

ChemicalEnvironmentConfig load_chemical_environment_config(const std::string &path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error(
            "Cannot open chemical environment config: " + path);

    json root;
    try
    {
        root = json::parse(in, /*callback=*/nullptr, /*allow_exceptions=*/true,
                           /*ignore_comments=*/true);
    }
    catch (const json::parse_error &e)
    {
        throw std::runtime_error(std::string("Invalid JSON in ") + path + ": " +
                                 e.what());
    }

    if (!root.contains("chemistry"))
        throw std::runtime_error(
            "chemical environment config: missing 'chemistry' section in " +
            path);

    const json &section = root.at("chemistry");

    ChemicalEnvironmentConfig cfg;

    if (section.contains("model"))
        cfg.model = section.at("model").get<std::string>();

    if (!section.contains("tick_interval_minutes"))
        throw std::runtime_error(
            "chemical environment config: tick_interval_minutes required");
    cfg.tick_interval_minutes = section.at("tick_interval_minutes").get<double>();
    if (cfg.tick_interval_minutes <= 0.0)
        throw std::invalid_argument(
            "chemical environment config: tick_interval_minutes must be > 0");

    if (!section.contains("channels"))
        throw std::runtime_error("chemical environment config: channels required");
    const json &channels = section.at("channels");
    cfg.channel_count = channels.at("count").get<int>();
    cfg.chemotaxis_channel = channels.at("chemotaxis").get<int>();
    if (cfg.channel_count <= 0)
        throw std::invalid_argument(
            "chemical environment config: channels.count must be positive");
    if (cfg.chemotaxis_channel < 0 ||
        cfg.chemotaxis_channel >= cfg.channel_count)
        throw std::invalid_argument(
            "chemical environment config: invalid channels.chemotaxis index");

    if (!section.contains("merge"))
        throw std::runtime_error("chemical environment config: merge required");
    cfg.merge_chemotaxis_from_species =
        section.at("merge").at("chemotaxis_from_species").get<std::string>();

    if (section.contains("baseline_total_mass"))
    {
        for (auto it = section.at("baseline_total_mass").begin();
             it != section.at("baseline_total_mass").end(); ++it)
            cfg.baseline_total_mass[it.key()] = it.value().get<float>();
    }

    if (!section.contains("species") || !section.at("species").is_array())
        throw std::runtime_error("chemical environment config: species array required");

    for (const json &entry : section.at("species"))
    {
        SpeciesConfigEntry s;
        s.id = entry.at("id").get<SpeciesId>();
        s.name = entry.at("name").get<std::string>();
        s.base_diffusivity_mm2_per_min =
            entry.at("base_diffusivity_mm2_per_min").get<double>();
        s.concentration_channel = entry.at("concentration_channel").get<int>();
        s.diffused_channel = entry.at("diffused_channel").get<int>();

        if (s.base_diffusivity_mm2_per_min <= 0.0)
            throw std::invalid_argument(
                "chemical environment config: base_diffusivity must be > 0 for " +
                s.name);
        if (s.concentration_channel < 0 || s.diffused_channel < 0)
            throw std::invalid_argument(
                "chemical environment config: channel indices must be >= 0 for " +
                s.name);
        if (s.concentration_channel >= cfg.channel_count ||
            s.diffused_channel >= cfg.channel_count)
            throw std::invalid_argument(
                "chemical environment config: channel index out of range for " +
                s.name);

        s.diffusivity_model = parse_diffusivity_model(entry, s.name);

        cfg.species.push_back(s);
    }

    for (const SpeciesConfigEntry &s : cfg.species)
        cfg.baseline_total_mass_for(s.name);

    if (cfg.species.empty())
        throw std::runtime_error(
            "chemical environment config: at least one species required");

    return cfg;
}
