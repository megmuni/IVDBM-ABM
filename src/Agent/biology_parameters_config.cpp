#include "biology_parameters_config.h"

#include "../World/Usr_World/biomaterialWorld.h"
#include "Usr_Agents/Cell.h"

#include <nlohmann/json.hpp>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

using json = nlohmann::json;

void reject_unknown_keys(const json &actual, const json &allowed,
                         const std::string &context) {
  for (auto it = actual.begin(); it != actual.end(); ++it) {
    if (it.key() == "description")
      continue;
    if (!allowed.contains(it.key()))
      throw std::invalid_argument("biology parameters config: unknown key '" +
                                  it.key() + "' in " + context);
    if (it.value().is_object() && allowed.at(it.key()).is_object())
      reject_unknown_keys(it.value(), allowed.at(it.key()),
                          context + "." + it.key());
  }
}

template <typename T>
T parse_strict(const json &obj, const std::string &context) {
  reject_unknown_keys(obj, json(T{}), context);
  return obj.get<T>();
}

void log_json_leaves(const json &node, const std::string &prefix) {
  for (auto it = node.begin(); it != node.end(); ++it) {
    const std::string key = prefix.empty() ? it.key() : prefix + "." + it.key();
    if (it.value().is_object())
      log_json_leaves(it.value(), key);
    else
      std::cout << "[biology_parameters] " << key << " = " << it.value()
                << std::endl;
  }
}

void flatten_json_leaves(const json &node, const std::string &prefix,
                         json &out) {
  for (auto it = node.begin(); it != node.end(); ++it) {
    const std::string key = prefix.empty() ? it.key() : prefix + "." + it.key();
    if (it.value().is_object())
      flatten_json_leaves(it.value(), key, out);
    else
      out[key] = it.value();
  }
}

} // namespace

BiologyParametersConfig
load_biology_parameters_config(const std::string &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open biology parameters config: " + path);

  json root;
  try {
    root = json::parse(in, /*callback=*/nullptr, /*allow_exceptions=*/true,
                       /*ignore_comments=*/true);
  } catch (const json::parse_error &e) {
    throw std::runtime_error(std::string("Invalid JSON in ") + path + ": " +
                             e.what());
  }

  if (!root.contains("biology"))
    throw std::runtime_error(
        "biology parameters config: missing 'biology' section in " + path);

  return parse_strict<BiologyParametersConfig>(root.at("biology"), "biology");
}

void apply_biology_parameters(const BiologyParametersConfig &cfg) {
  const auto &cp = cfg.cell.proliferation;
  Cell::proliferation[Cell::PROLIFERATION_HOURS_BETWEEN] = static_cast<float>(cp.hours_between_proliferation);
  Cell::proliferation[Cell::PROLIFERATION_TGF_THRESHOLD] = static_cast<float>(cp.tgf_threshold);
  Cell::proliferation[Cell::PROLIFERATION_LOG_SCALE] = static_cast<float>(cp.log_scale);
  Cell::proliferation[Cell::PROLIFERATION_LOG_OFFSET] = static_cast<float>(cp.log_offset);

  const auto &cs = cfg.cell.cytokine_synthesis;
  Cell::cytokineSynthesis[Cell::CYTOKINE_TGF_BASELINE] = static_cast<float>(cs.tgf_baseline);
  Cell::cytokineSynthesis[Cell::CYTOKINE_TGF_FEEDBACK_TGF] = static_cast<float>(cs.tgf_feedback_tgf);
  Cell::cytokineSynthesis[Cell::CYTOKINE_TGF_FEEDBACK_IL1BETA] = static_cast<float>(cs.tgf_feedback_il1beta);
  Cell::cytokineSynthesis[Cell::CYTOKINE_TGF_FEEDBACK_TNF] = static_cast<float>(cs.tgf_feedback_tnf);
  Cell::cytokineSynthesis[Cell::CYTOKINE_TNF_BASELINE] = static_cast<float>(cs.tnf_baseline);
  Cell::cytokineSynthesis[Cell::CYTOKINE_TNF_FEEDBACK_IL1BETA] = static_cast<float>(cs.tnf_feedback_il1beta);
  Cell::cytokineSynthesis[Cell::CYTOKINE_TNF_FEEDBACK_TGF_DENOM] = static_cast<float>(cs.tnf_feedback_tgf_denom);
  Cell::cytokineSynthesis[Cell::CYTOKINE_IL1BETA_BASELINE] = static_cast<float>(cs.il1beta_baseline);
  Cell::cytokineSynthesis[Cell::CYTOKINE_IL1BETA_FEEDBACK_TNF] = static_cast<float>(cs.il1beta_feedback_tnf);
  Cell::cytokineSynthesis[Cell::CYTOKINE_IL1BETA_FEEDBACK_TGF_DENOM] =
      static_cast<float>(cs.il1beta_feedback_tgf_denom);

  const auto &st = cfg.stem;
  Stem::OCR = static_cast<float>(st.ocr_fmol_per_hour_per_cell / 2.0);
  Stem::apoptosisChance = static_cast<float>(st.apoptosis_chance);
  Stem::CaAlgMigration[Stem::MIGRATION_ELASTICITY_EFFECT] = static_cast<float>(st.migration.elasticity_effect);
  Stem::CaAlgMigration[Stem::MIGRATION_BASELINE_SPEED] = static_cast<float>(st.migration.baseline_speed);
  Stem::cytokineSynthesis[Stem::CYTOKINE_TGF_BASELINE] =
      static_cast<float>(st.cytokine_synthesis.tgf_baseline);
  Stem::cytokineSynthesis[Stem::CYTOKINE_TNF_BASELINE] =
      static_cast<float>(st.cytokine_synthesis.tnf_baseline);
  Stem::cytokineSynthesis[Stem::CYTOKINE_IL1BETA_BASELINE] =
      static_cast<float>(st.cytokine_synthesis.il1beta_baseline);
  Stem::CollagenSynth[0] =
      static_cast<float>(st.collagen_synthesis.baseline_rate);
  Stem::AggrecanSynth[0] =
      static_cast<float>(st.aggrecan_synthesis.tgf_threshold);
  Stem::proliferation[Stem::PROLIFERATION_TGF_THRESHOLD] = static_cast<float>(st.proliferation.tgf_threshold);
  Stem::proliferation[Stem::PROLIFERATION_TNF_EFFECT] = static_cast<float>(st.proliferation.tnf_effect);
  Stem::proliferation[Stem::PROLIFERATION_IL1BETA_EFFECT] = static_cast<float>(st.proliferation.il1beta_effect);
  Stem::proliferation[Stem::PROLIFERATION_ELASTICITY_EFFECT] =
      static_cast<float>(st.proliferation.elasticity_effect);
  Stem::differentiation[Stem::DIFFERENTIATION_ASYMMETRIC_PROBABILITY] =
      static_cast<float>(st.differentiation.asymmetric_probability);
  Stem::differentiation[Stem::DIFFERENTIATION_BASELINE_PROBABILITY] =
      static_cast<float>(st.differentiation.baseline_probability);
  Stem::differentiation[Stem::DIFFERENTIATION_TGF_EFFECT] = static_cast<float>(st.differentiation.tgf_effect);
  Stem::differentiation[Stem::DIFFERENTIATION_HOURS_BETWEEN_ATTEMPTS] =
      static_cast<float>(st.differentiation.hours_between_attempts);

  const auto &pg = cfg.progen;
  Progen::OCR = static_cast<float>(pg.ocr_fmol_per_hour_per_cell / 2.0);
  Progen::apoptosisChance = static_cast<float>(pg.apoptosis_chance);
  Progen::CaAlgMigration[Progen::MIGRATION_ELASTICITY_EFFECT] =
      static_cast<float>(pg.migration.elasticity_effect);
  Progen::CaAlgMigration[Progen::MIGRATION_BASELINE_SPEED] = static_cast<float>(pg.migration.baseline_speed);
  Progen::cytokineSynthesis[Progen::CYTOKINE_TGF_BASELINE] =
      static_cast<float>(pg.cytokine_synthesis.tgf_baseline);
  Progen::cytokineSynthesis[Progen::CYTOKINE_TNF_BASELINE] =
      static_cast<float>(pg.cytokine_synthesis.tnf_baseline);
  Progen::cytokineSynthesis[Progen::CYTOKINE_IL1BETA_BASELINE] =
      static_cast<float>(pg.cytokine_synthesis.il1beta_baseline);
  Progen::AggrecanSynth[0] =
      static_cast<float>(pg.aggrecan_synthesis.baseline_rate);

  const auto &np = cfg.np;
  NP::OCR = static_cast<float>(np.ocr_fmol_per_hour_per_cell / 2.0);
  NP::apoptosisChance = static_cast<float>(np.apoptosis_chance);
  NP::CaAlgMigration[NP::MIGRATION_ELASTICITY_EFFECT] = static_cast<float>(np.migration.elasticity_effect);
  NP::CaAlgMigration[NP::MIGRATION_BASELINE_SPEED] = static_cast<float>(np.migration.baseline_speed);
  NP::CollagenSynth[NP::COLLAGEN_SCALING_FACTOR] =
      static_cast<float>(np.collagen_synthesis.scaling_factor);
  NP::CollagenSynth[NP::COLLAGEN_TIME_EFFECT] = static_cast<float>(np.collagen_synthesis.time_effect);
  NP::CollagenSynth[NP::COLLAGEN_BASELINE_RATE] =
      static_cast<float>(np.collagen_synthesis.baseline_rate);
  NP::AggrecanSynth[NP::AGGRECAN_SCALING_FACTOR] =
      static_cast<float>(np.aggrecan_synthesis.scaling_factor);
  NP::AggrecanSynth[NP::AGGRECAN_TIME_EFFECT] = static_cast<float>(np.aggrecan_synthesis.time_effect);
  NP::AggrecanSynth[NP::AGGRECAN_BASELINE_RATE] =
      static_cast<float>(np.aggrecan_synthesis.baseline_rate);

  const auto &bm = cfg.biomaterial;
  BMWorld::ElasticMod[BMWorld::ELASTIC_INTERCEPT] = static_cast<float>(bm.elastic_modulus.intercept);
  BMWorld::ElasticMod[BMWorld::ELASTIC_ALGINATE_CONCENTRATION] =
      static_cast<float>(bm.elastic_modulus.alginate_concentration);
  BMWorld::ElasticMod[BMWorld::ELASTIC_CROSSLINKER_DENSITY] =
      static_cast<float>(bm.elastic_modulus.crosslinker_density);
  BMWorld::ElasticMod[BMWorld::ELASTIC_ALGINATE_MOLECULAR_WEIGHT] =
      static_cast<float>(bm.elastic_modulus.alginate_molecular_weight);
  BMWorld::ElasticMod[BMWorld::ELASTIC_ALGINATE_CROSSLINKER_INTERACTION] =
      static_cast<float>(bm.elastic_modulus.alginate_crosslinker_interaction);
  BMWorld::ElasticMod[BMWorld::ELASTIC_ALGINATE_MW_INTERACTION] =
      static_cast<float>(bm.elastic_modulus.alginate_mw_interaction);
  BMWorld::ElasticMod[BMWorld::ELASTIC_MW_CROSSLINKER_INTERACTION] =
      static_cast<float>(bm.elastic_modulus.mw_crosslinker_interaction);
  BMWorld::PoreSize[BMWorld::PORE_CROSSLINKER_EFFECT] = static_cast<float>(bm.pore_size.crosslinker_effect);
  BMWorld::PoreSize[BMWorld::PORE_BASELINE] = static_cast<float>(bm.pore_size.baseline);
  BMWorld::MassLoss[BMWorld::MASSLOSS_BASELINE] = static_cast<float>(bm.mass_loss.baseline);
  BMWorld::MassLoss[BMWorld::MASSLOSS_CROSSLINKER_EFFECT] = static_cast<float>(bm.mass_loss.crosslinker_effect);
  BMWorld::MassLoss[BMWorld::MASSLOSS_TIME_EFFECT] = static_cast<float>(bm.mass_loss.time_effect);
  BMWorld::MassLoss[BMWorld::MASSLOSS_CROSSLINKER_TIME_INTERACTION] =
      static_cast<float>(bm.mass_loss.crosslinker_time_interaction);
  BMWorld::SwellRatio[BMWorld::SWELL_BASELINE] = static_cast<float>(bm.swell_ratio.baseline);
  BMWorld::SwellRatio[BMWorld::SWELL_TIME_EFFECT] = static_cast<float>(bm.swell_ratio.time_effect);
  BMWorld::SwellRatio[BMWorld::SWELL_ALGINATE_CONCENTRATION_EFFECT] =
      static_cast<float>(bm.swell_ratio.alginate_concentration_effect);
  BMWorld::SwellRatio[BMWorld::SWELL_TIME_CROSSLINKER_INTERACTION] =
      static_cast<float>(bm.swell_ratio.time_crosslinker_interaction);
  BMWorld::SwellRatio[BMWorld::SWELL_ALGINATE_CROSSLINKER_INTERACTION] =
      static_cast<float>(bm.swell_ratio.alginate_crosslinker_interaction);
}

void log_biology_parameters(const BiologyParametersConfig &cfg,
                            const std::string &source_path) {
  std::cout << "[biology_parameters] source: " << source_path << std::endl;
  log_json_leaves(json(cfg), "");
}

void record_biology_parameters_in_run_params(const BiologyParametersConfig &cfg,
                                             const std::string &source_path) {
  const char *path = std::getenv("IVDBM_RUN_PARAMS_JSON");
  if (path == nullptr || path[0] == '\0')
    return;

  json root;
  {
    std::ifstream in(path);
    if (in)
      in >> root;
  }

  json values = json::object();
  flatten_json_leaves(json(cfg), "", values);

  if (!root.contains("simulation"))
    root["simulation"] = json::object();
  root["simulation"]["simulation_config"] = source_path;
  root["biology_parameters"] = values;

  std::ofstream out(path);
  if (!out) {
    std::fprintf(stderr, "Warning: cannot update run params at %s\n", path);
    return;
  }
  out << root.dump(2) << std::endl;
}
