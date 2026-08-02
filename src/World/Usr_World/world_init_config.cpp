#include "world_init_config.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <set>
#include <stdexcept>

namespace
{

using json = nlohmann::json;

void reject_unknown_keys(const json &obj, const std::set<std::string> &allowed,
                         const std::string &context)
{
  for (auto it = obj.begin(); it != obj.end(); ++it)
  {
    if (allowed.count(it.key()) == 0)
      throw std::invalid_argument("world_init config: unknown key '" +
                                  it.key() + "' in " + context);
  }
}

} // namespace

WorldInitParams load_world_init_config(const std::string &path)
{
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open simulation config: " + path);

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

  if (!root.contains("world_init"))
    throw std::runtime_error(
        "world_init config: missing 'world_init' section in " + path);

  const json &section = root.at("world_init");
  reject_unknown_keys(section,
                      {"description", "msc_count", "alginate", "peptide"},
                      "world_init");

  WorldInitParams cfg;
  cfg.msc_count = section.at("msc_count").get<int>();

  const json &alg = section.at("alginate");
  reject_unknown_keys(
      alg, {"wv_percent", "high_mw_ratio", "low_mw_ratio", "ca_mm"},
      "world_init.alginate");
  cfg.alginate.wv_percent = alg.at("wv_percent").get<double>();
  cfg.alginate.high_mw_ratio = alg.at("high_mw_ratio").get<double>();
  cfg.alginate.low_mw_ratio = alg.at("low_mw_ratio").get<double>();
  cfg.alginate.ca_mm = alg.at("ca_mm").get<double>();

  if (section.contains("peptide"))
    cfg.peptide = section.at("peptide").get<std::string>();

  return cfg;
}
