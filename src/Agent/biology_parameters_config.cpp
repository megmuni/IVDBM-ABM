#include "biology_parameters_config.h"

#include "Usr_Agents/Cell.h"
#include "../World/Usr_World/biomaterialWorld.h"

#include <nlohmann/json.hpp>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
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
      throw std::invalid_argument("biology parameters config: unknown key '" +
                                  it.key() + "' in " + context);
  }
}

double require_double(const json &obj, const char *key, const std::string &context)
{
  if (!obj.contains(key))
    throw std::invalid_argument("biology parameters config: missing '" +
                                std::string(key) + "' in " + context);
  return obj.at(key).get<double>();
}

MigrationParams parse_migration(const json &obj, const std::string &context)
{
  reject_unknown_keys(obj, {"description", "elasticity_effect", "baseline_speed"},
                      context);
  MigrationParams m;
  m.elasticity_effect = require_double(obj, "elasticity_effect", context);
  m.baseline_speed = require_double(obj, "baseline_speed", context);
  return m;
}

CellProliferationParams parse_cell_proliferation(const json &obj)
{
  reject_unknown_keys(obj,
                      {"description", "hours_between_proliferation", "tgf_threshold",
                       "log_scale", "log_offset"},
                      "cell.proliferation");
  CellProliferationParams p;
  p.hours_between_proliferation =
      require_double(obj, "hours_between_proliferation", "cell.proliferation");
  p.tgf_threshold = require_double(obj, "tgf_threshold", "cell.proliferation");
  p.log_scale = require_double(obj, "log_scale", "cell.proliferation");
  p.log_offset = require_double(obj, "log_offset", "cell.proliferation");
  return p;
}

CellCytokineSynthesisParams parse_cell_cytokine_synthesis(const json &obj)
{
  reject_unknown_keys(
      obj,
      {"description", "tgf_baseline", "tgf_feedback_tgf", "tgf_feedback_il1beta",
       "tgf_feedback_tnf", "tnf_baseline", "tnf_feedback_il1beta",
       "tnf_feedback_tgf_denom", "il1beta_baseline", "il1beta_feedback_tnf",
       "il1beta_feedback_tgf_denom"},
      "cell.cytokine_synthesis");
  CellCytokineSynthesisParams c;
  c.tgf_baseline =
      require_double(obj, "tgf_baseline", "cell.cytokine_synthesis");
  c.tgf_feedback_tgf =
      require_double(obj, "tgf_feedback_tgf", "cell.cytokine_synthesis");
  c.tgf_feedback_il1beta =
      require_double(obj, "tgf_feedback_il1beta", "cell.cytokine_synthesis");
  c.tgf_feedback_tnf =
      require_double(obj, "tgf_feedback_tnf", "cell.cytokine_synthesis");
  c.tnf_baseline =
      require_double(obj, "tgf_baseline", "cell.cytokine_synthesis");
  c.tnf_feedback_il1beta =
      require_double(obj, "tnf_feedback_il1beta", "cell.cytokine_synthesis");
  c.tnf_feedback_tgf_denom =
      require_double(obj, "tnf_feedback_tgf_denom", "cell.cytokine_synthesis");
  c.il1beta_baseline =
      require_double(obj, "il1beta_baseline", "cell.cytokine_synthesis");
  c.il1beta_feedback_tnf =
      require_double(obj, "il1beta_feedback_tnf", "cell.cytokine_synthesis");
  c.il1beta_feedback_tgf_denom =
      require_double(obj, "il1beta_feedback_tgf_denom", "cell.cytokine_synthesis");
  return c;
}

StemParams parse_stem(const json &obj)
{
  reject_unknown_keys(
      obj,
      {"ocr_fmol_per_hour_per_cell", "apoptosis_chance", "migration",
       "cytokine_synthesis", "collagen_synthesis", "aggrecan_synthesis",
       "proliferation", "differentiation"},
      "stem");
  StemParams s;
  s.ocr_fmol_per_hour_per_cell =
      require_double(obj, "ocr_fmol_per_hour_per_cell", "stem");
  s.apoptosis_chance = require_double(obj, "apoptosis_chance", "stem");
  s.migration = parse_migration(obj.at("migration"), "stem.migration");

  const json &cyto = obj.at("cytokine_synthesis");
  reject_unknown_keys(cyto, {"description", "tgf_baseline", "tnf_baseline",
                             "il1beta_baseline"},
                      "stem.cytokine_synthesis");
  s.cytokine_synthesis.tgf_baseline =
      require_double(cyto, "tgf_baseline", "stem.cytokine_synthesis");
  s.cytokine_synthesis.tnf_baseline =
      require_double(cyto, "tnf_baseline", "stem.cytokine_synthesis");
  s.cytokine_synthesis.il1beta_baseline =
      require_double(cyto, "il1beta_baseline", "stem.cytokine_synthesis");

  const json &coll = obj.at("collagen_synthesis");
  reject_unknown_keys(coll, {"baseline_rate"}, "stem.collagen_synthesis");
  s.collagen_synthesis_baseline_rate =
      require_double(coll, "baseline_rate", "stem.collagen_synthesis");

  const json &agg = obj.at("aggrecan_synthesis");
  reject_unknown_keys(agg, {"tgf_threshold"}, "stem.aggrecan_synthesis");
  s.aggrecan_synthesis_tgf_threshold =
      require_double(agg, "tgf_threshold", "stem.aggrecan_synthesis");

  const json &prolif = obj.at("proliferation");
  reject_unknown_keys(prolif,
                      {"tgf_threshold", "tnf_effect", "il1beta_effect",
                       "elasticity_effect"},
                      "stem.proliferation");
  s.proliferation.tgf_threshold =
      require_double(prolif, "tgf_threshold", "stem.proliferation");
  s.proliferation.tnf_effect =
      require_double(prolif, "tnf_effect", "stem.proliferation");
  s.proliferation.il1beta_effect =
      require_double(prolif, "il1beta_effect", "stem.proliferation");
  s.proliferation.elasticity_effect =
      require_double(prolif, "elasticity_effect", "stem.proliferation");

  const json &diff = obj.at("differentiation");
  reject_unknown_keys(
      diff,
      {"asymmetric_probability", "baseline_probability", "tgf_effect",
       "hours_between_attempts"},
      "stem.differentiation");
  s.differentiation.asymmetric_probability =
      require_double(diff, "asymmetric_probability", "stem.differentiation");
  s.differentiation.baseline_probability =
      require_double(diff, "baseline_probability", "stem.differentiation");
  s.differentiation.tgf_effect =
      require_double(diff, "tgf_effect", "stem.differentiation");
  s.differentiation.hours_between_attempts =
      require_double(diff, "hours_between_attempts", "stem.differentiation");
  return s;
}

ProgenParams parse_progen(const json &obj)
{
  reject_unknown_keys(obj,
                      {"ocr_fmol_per_hour_per_cell", "apoptosis_chance",
                       "migration", "cytokine_synthesis", "aggrecan_synthesis"},
                      "progen");
  ProgenParams p;
  p.ocr_fmol_per_hour_per_cell =
      require_double(obj, "ocr_fmol_per_hour_per_cell", "progen");
  p.apoptosis_chance = require_double(obj, "apoptosis_chance", "progen");
  p.migration = parse_migration(obj.at("migration"), "progen.migration");

  const json &cyto = obj.at("cytokine_synthesis");
  reject_unknown_keys(cyto, {"tgf_baseline", "tnf_baseline", "il1beta_baseline"},
                      "progen.cytokine_synthesis");
  p.cytokine_synthesis.tgf_baseline =
      require_double(cyto, "tgf_baseline", "progen.cytokine_synthesis");
  p.cytokine_synthesis.tnf_baseline =
      require_double(cyto, "tnf_baseline", "progen.cytokine_synthesis");
  p.cytokine_synthesis.il1beta_baseline =
      require_double(cyto, "il1beta_baseline", "progen.cytokine_synthesis");

  const json &agg = obj.at("aggrecan_synthesis");
  reject_unknown_keys(agg, {"baseline_rate"}, "progen.aggrecan_synthesis");
  p.aggrecan_synthesis_baseline_rate =
      require_double(agg, "baseline_rate", "progen.aggrecan_synthesis");
  return p;
}

NpParams parse_np(const json &obj)
{
  reject_unknown_keys(obj,
                      {"ocr_fmol_per_hour_per_cell", "apoptosis_chance",
                       "migration", "collagen_synthesis", "aggrecan_synthesis"},
                      "np");
  NpParams n;
  n.ocr_fmol_per_hour_per_cell =
      require_double(obj, "ocr_fmol_per_hour_per_cell", "np");
  n.apoptosis_chance = require_double(obj, "apoptosis_chance", "np");
  n.migration = parse_migration(obj.at("migration"), "np.migration");

  const json &coll = obj.at("collagen_synthesis");
  reject_unknown_keys(coll, {"scaling_factor", "time_effect", "baseline_rate"},
                      "np.collagen_synthesis");
  n.collagen_synthesis.scaling_factor =
      require_double(coll, "scaling_factor", "np.collagen_synthesis");
  n.collagen_synthesis.time_effect =
      require_double(coll, "time_effect", "np.collagen_synthesis");
  n.collagen_synthesis.baseline_rate =
      require_double(coll, "baseline_rate", "np.collagen_synthesis");

  const json &agg = obj.at("aggrecan_synthesis");
  reject_unknown_keys(agg, {"scaling_factor", "time_effect", "baseline_rate"},
                      "np.aggrecan_synthesis");
  n.aggrecan_synthesis.scaling_factor =
      require_double(agg, "scaling_factor", "np.aggrecan_synthesis");
  n.aggrecan_synthesis.time_effect =
      require_double(agg, "time_effect", "np.aggrecan_synthesis");
  n.aggrecan_synthesis.baseline_rate =
      require_double(agg, "baseline_rate", "np.aggrecan_synthesis");
  return n;
}

BiomaterialParams parse_biomaterial(const json &obj)
{
  reject_unknown_keys(obj, {"elastic_modulus", "pore_size", "mass_loss",
                            "swell_ratio"},
                      "biomaterial");
  BiomaterialParams b;

  const json &em = obj.at("elastic_modulus");
  reject_unknown_keys(
      em,
      {"description", "intercept", "alginate_concentration", "crosslinker_density",
       "alginate_molecular_weight", "alginate_crosslinker_interaction",
       "alginate_mw_interaction", "mw_crosslinker_interaction"},
      "biomaterial.elastic_modulus");
  b.elastic_modulus.intercept =
      require_double(em, "intercept", "biomaterial.elastic_modulus");
  b.elastic_modulus.alginate_concentration =
      require_double(em, "alginate_concentration", "biomaterial.elastic_modulus");
  b.elastic_modulus.crosslinker_density =
      require_double(em, "crosslinker_density", "biomaterial.elastic_modulus");
  b.elastic_modulus.alginate_molecular_weight =
      require_double(em, "alginate_molecular_weight", "biomaterial.elastic_modulus");
  b.elastic_modulus.alginate_crosslinker_interaction =
      require_double(em, "alginate_crosslinker_interaction",
                       "biomaterial.elastic_modulus");
  b.elastic_modulus.alginate_mw_interaction =
      require_double(em, "alginate_mw_interaction", "biomaterial.elastic_modulus");
  b.elastic_modulus.mw_crosslinker_interaction =
      require_double(em, "mw_crosslinker_interaction", "biomaterial.elastic_modulus");

  const json &ps = obj.at("pore_size");
  reject_unknown_keys(ps, {"crosslinker_effect", "baseline"},
                      "biomaterial.pore_size");
  b.pore_size.crosslinker_effect =
      require_double(ps, "crosslinker_effect", "biomaterial.pore_size");
  b.pore_size.baseline = require_double(ps, "baseline", "biomaterial.pore_size");

  const json &ml = obj.at("mass_loss");
  reject_unknown_keys(ml,
                      {"baseline", "crosslinker_effect", "time_effect",
                       "crosslinker_time_interaction"},
                      "biomaterial.mass_loss");
  b.mass_loss.baseline = require_double(ml, "baseline", "biomaterial.mass_loss");
  b.mass_loss.crosslinker_effect =
      require_double(ml, "crosslinker_effect", "biomaterial.mass_loss");
  b.mass_loss.time_effect =
      require_double(ml, "time_effect", "biomaterial.mass_loss");
  b.mass_loss.crosslinker_time_interaction =
      require_double(ml, "crosslinker_time_interaction", "biomaterial.mass_loss");

  const json &sr = obj.at("swell_ratio");
  reject_unknown_keys(
      sr,
      {"baseline", "time_effect", "alginate_concentration_effect",
       "time_crosslinker_interaction", "alginate_crosslinker_interaction"},
      "biomaterial.swell_ratio");
  b.swell_ratio.baseline =
      require_double(sr, "baseline", "biomaterial.swell_ratio");
  b.swell_ratio.time_effect =
      require_double(sr, "time_effect", "biomaterial.swell_ratio");
  b.swell_ratio.alginate_concentration_effect =
      require_double(sr, "alginate_concentration_effect", "biomaterial.swell_ratio");
  b.swell_ratio.time_crosslinker_interaction =
      require_double(sr, "time_crosslinker_interaction", "biomaterial.swell_ratio");
  b.swell_ratio.alginate_crosslinker_interaction =
      require_double(sr, "alginate_crosslinker_interaction", "biomaterial.swell_ratio");
  return b;
}

void log_entry(const std::string &key, double value)
{
  std::cout << "[biology_parameters] " << key << " = " << value << std::endl;
}

} // namespace

BiologyParametersConfig load_biology_parameters_config(const std::string &path)
{
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open biology parameters config: " + path);

  json root;
  try
  {
    // ignore_comments=true so the template's "// previous label" annotations
    // (e.g. "m5", "c10") can stay next to each field without failing to parse.
    root = json::parse(in, /*callback=*/nullptr, /*allow_exceptions=*/true,
                       /*ignore_comments=*/true);
  }
  catch (const json::parse_error &e)
  {
    throw std::runtime_error(std::string("Invalid JSON in ") + path + ": " +
                             e.what());
  }

  if (!root.contains("biology"))
    throw std::runtime_error(
        "biology parameters config: missing 'biology' section in " + path);

  const json &section = root.at("biology");

  reject_unknown_keys(section,
                      {"description", "cell", "stem", "progen",
                       "np", "biomaterial"},
                      "biology");

  BiologyParametersConfig cfg;

  const json &cell = section.at("cell");
  reject_unknown_keys(cell, {"proliferation", "cytokine_synthesis"}, "cell");
  cfg.cell.proliferation = parse_cell_proliferation(cell.at("proliferation"));
  cfg.cell.cytokine_synthesis =
      parse_cell_cytokine_synthesis(cell.at("cytokine_synthesis"));

  cfg.stem = parse_stem(section.at("stem"));
  cfg.progen = parse_progen(section.at("progen"));
  cfg.np = parse_np(section.at("np"));
  cfg.biomaterial = parse_biomaterial(section.at("biomaterial"));

  return cfg;
}

void apply_biology_parameters(const BiologyParametersConfig &cfg)
{
  const auto &cp = cfg.cell.proliferation;
  // Hours between proliferation attempts; Cell::proliferate (Cell.cpp)
  Cell::proliferation[0] = static_cast<float>(cp.hours_between_proliferation);
  // TGF threshold for NP proliferation probability; NP::get_prolif_prob
  Cell::proliferation[1] = static_cast<float>(cp.tgf_threshold);
  // Log-scaling factor for NP proliferation probability; NP::get_prolif_prob
  Cell::proliferation[2] = static_cast<float>(cp.log_scale);
  // Offset for NP proliferation probability; NP::get_prolif_prob
  Cell::proliferation[3] = static_cast<float>(cp.log_offset);

  const auto &cs = cfg.cell.cytokine_synthesis;
  // Baseline TGF synthesis (shared feedback model); create_cytokines (Cell.cpp)
  Cell::cytokineSynthesis[0] = static_cast<float>(cs.tgf_baseline);
  // Effect of local TGF on TGF synthesis; create_cytokines
  Cell::cytokineSynthesis[1] = static_cast<float>(cs.tgf_feedback_tgf);
  // Effect of IL-1β on TGF synthesis; create_cytokines
  Cell::cytokineSynthesis[2] = static_cast<float>(cs.tgf_feedback_il1beta);
  // Effect of TNF on TGF synthesis; create_cytokines
  Cell::cytokineSynthesis[3] = static_cast<float>(cs.tgf_feedback_tnf);
  // Baseline TNF synthesis; create_cytokines
  Cell::cytokineSynthesis[4] = static_cast<float>(cs.tnf_baseline);
  // Effect of IL-1β on TNF synthesis; create_cytokines
  Cell::cytokineSynthesis[5] = static_cast<float>(cs.tnf_feedback_il1beta);
  // TGF denominator in TNF synthesis feedback; create_cytokines
  Cell::cytokineSynthesis[6] = static_cast<float>(cs.tnf_feedback_tgf_denom);
  // Baseline IL-1β synthesis; create_cytokines
  Cell::cytokineSynthesis[7] = static_cast<float>(cs.il1beta_baseline);
  // Effect of TNF on IL-1β synthesis; create_cytokines
  Cell::cytokineSynthesis[8] = static_cast<float>(cs.il1beta_feedback_tnf);
  // TGF denominator in IL-1β synthesis feedback; create_cytokines
  Cell::cytokineSynthesis[9] = static_cast<float>(cs.il1beta_feedback_tgf_denom);

  const auto &st = cfg.stem;
  // MSC oxygen consumption (fmol/h/cell → per-tick); Stem::get_OCR
  Stem::OCR = static_cast<float>(st.ocr_fmol_per_hour_per_cell / 2.0);
  // MSC apoptosis probability; Stem::get_apoptosis_chance
  Stem::apoptosisChance = static_cast<float>(st.apoptosis_chance);
  // Effect of hydrogel elasticity on MSC migration; Stem::get_migration_speed
  Stem::CaAlgMigration[0] = static_cast<float>(st.migration.elasticity_effect);
  // Baseline MSC migration speed; Stem::get_migration_speed
  Stem::CaAlgMigration[1] = static_cast<float>(st.migration.baseline_speed);
  // Baseline MSC TGF secretion; Stem::create_cytokines
  Stem::cytokineSynthesis[0] = static_cast<float>(st.cytokine_synthesis.tgf_baseline);
  // Baseline MSC TNF secretion; Stem::create_cytokines
  Stem::cytokineSynthesis[1] = static_cast<float>(st.cytokine_synthesis.tnf_baseline);
  // Baseline MSC IL-1β secretion; Stem::create_cytokines
  Stem::cytokineSynthesis[2] =
      static_cast<float>(st.cytokine_synthesis.il1beta_baseline);
  // Baseline MSC collagen synthesis rate; Stem::calculate_ecm_synth_rates
  Stem::CollagenSynth[0] = static_cast<float>(st.collagen_synthesis_baseline_rate);
  // TGF threshold for MSC aggrecan synthesis; Stem::calculate_ecm_synth_rates
  Stem::AggrecanSynth[0] =
      static_cast<float>(st.aggrecan_synthesis_tgf_threshold);
  // TGF threshold for MSC proliferation; Stem::get_prolif_prob
  Stem::proliferation[0] = static_cast<float>(st.proliferation.tgf_threshold);
  // Effect of TNF on MSC proliferation (commented formula); Stem::get_prolif_prob
  Stem::proliferation[1] = static_cast<float>(st.proliferation.tnf_effect);
  // Effect of IL-1β on MSC proliferation; Stem::get_prolif_prob
  Stem::proliferation[2] = static_cast<float>(st.proliferation.il1beta_effect);
  // Effect of elasticity on MSC proliferation; Stem::get_prolif_prob
  Stem::proliferation[3] = static_cast<float>(st.proliferation.elasticity_effect);
  // Probability of asymmetric MSC differentiation; Cell::differentiate
  Stem::differentiation[0] =
      static_cast<float>(st.differentiation.asymmetric_probability);
  // Baseline MSC differentiation probability; Stem::get_diff_prob
  Stem::differentiation[1] =
      static_cast<float>(st.differentiation.baseline_probability);
  // Effect of TGF on differentiation probability; Stem::get_diff_prob
  Stem::differentiation[2] = static_cast<float>(st.differentiation.tgf_effect);
  // Hours between MSC differentiation attempts; Cell::differentiate
  Stem::differentiation[3] =
      static_cast<float>(st.differentiation.hours_between_attempts);

  const auto &pg = cfg.progen;
  // Pre-NP oxygen consumption (fmol/h/cell → per-tick); Progen::get_OCR
  Progen::OCR = static_cast<float>(pg.ocr_fmol_per_hour_per_cell / 2.0);
  // Pre-NP apoptosis probability; Progen::get_apoptosis_chance
  Progen::apoptosisChance = static_cast<float>(pg.apoptosis_chance);
  // Effect of elasticity on Pre-NP migration; Progen::get_migration_speed
  Progen::CaAlgMigration[0] = static_cast<float>(pg.migration.elasticity_effect);
  // Baseline Pre-NP migration speed; Progen::get_migration_speed
  Progen::CaAlgMigration[1] = static_cast<float>(pg.migration.baseline_speed);
  // Baseline Pre-NP TGF secretion; Progen::create_cytokines
  Progen::cytokineSynthesis[0] =
      static_cast<float>(pg.cytokine_synthesis.tgf_baseline);
  // Baseline Pre-NP TNF secretion; Progen::create_cytokines
  Progen::cytokineSynthesis[1] =
      static_cast<float>(pg.cytokine_synthesis.tnf_baseline);
  // Baseline Pre-NP IL-1β secretion; Progen::create_cytokines
  Progen::cytokineSynthesis[2] =
      static_cast<float>(pg.cytokine_synthesis.il1beta_baseline);
  // Baseline Pre-NP aggrecan synthesis rate; Progen::calculate_ecm_synth_rates
  Progen::AggrecanSynth[0] =
      static_cast<float>(pg.aggrecan_synthesis_baseline_rate);

  const auto &np = cfg.np;
  // NP oxygen consumption (fmol/h/cell → per-tick); NP::get_OCR
  NP::OCR = static_cast<float>(np.ocr_fmol_per_hour_per_cell / 2.0);
  // NP apoptosis probability; NP::get_apoptosis_chance
  NP::apoptosisChance = static_cast<float>(np.apoptosis_chance);
  // Effect of elasticity on NP migration; NP::get_migration_speed
  NP::CaAlgMigration[0] = static_cast<float>(np.migration.elasticity_effect);
  // Baseline NP migration speed; NP::get_migration_speed
  NP::CaAlgMigration[1] = static_cast<float>(np.migration.baseline_speed);
  // NP collagen synthesis scaling factor; NP::calculate_ecm_synth_rates
  NP::CollagenSynth[0] =
      static_cast<float>(np.collagen_synthesis.scaling_factor);
  // Effect of culture time on NP collagen synthesis; NP::calculate_ecm_synth_rates
  NP::CollagenSynth[1] = static_cast<float>(np.collagen_synthesis.time_effect);
  // Baseline NP collagen synthesis rate; NP::calculate_ecm_synth_rates
  NP::CollagenSynth[2] =
      static_cast<float>(np.collagen_synthesis.baseline_rate);
  // NP aggrecan synthesis scaling factor; NP::calculate_ecm_synth_rates
  NP::AggrecanSynth[0] =
      static_cast<float>(np.aggrecan_synthesis.scaling_factor);
  // Effect of culture time on NP aggrecan synthesis; NP::calculate_ecm_synth_rates
  NP::AggrecanSynth[1] = static_cast<float>(np.aggrecan_synthesis.time_effect);
  // Baseline NP aggrecan synthesis rate; NP::calculate_ecm_synth_rates
  NP::AggrecanSynth[2] =
      static_cast<float>(np.aggrecan_synthesis.baseline_rate);

  const auto &bm = cfg.biomaterial;
  // Elastic modulus intercept; BMWorld::initializeCaAlg
  BMWorld::ElasticMod[0] = static_cast<float>(bm.elastic_modulus.intercept);
  // Effect of alginate concentration on E; BMWorld::initializeCaAlg
  BMWorld::ElasticMod[1] =
      static_cast<float>(bm.elastic_modulus.alginate_concentration);
  // Effect of crosslinker density on E; BMWorld::initializeCaAlg
  BMWorld::ElasticMod[2] =
      static_cast<float>(bm.elastic_modulus.crosslinker_density);
  // Effect of alginate molecular weight on E; BMWorld::initializeCaAlg
  BMWorld::ElasticMod[3] =
      static_cast<float>(bm.elastic_modulus.alginate_molecular_weight);
  // Alginate × crosslinker interaction on E; BMWorld::initializeCaAlg
  BMWorld::ElasticMod[4] =
      static_cast<float>(bm.elastic_modulus.alginate_crosslinker_interaction);
  // Alginate × molecular weight interaction on E; BMWorld::initializeCaAlg
  BMWorld::ElasticMod[5] =
      static_cast<float>(bm.elastic_modulus.alginate_mw_interaction);
  // Molecular weight × crosslinker interaction on E; BMWorld::initializeCaAlg
  BMWorld::ElasticMod[6] =
      static_cast<float>(bm.elastic_modulus.mw_crosslinker_interaction);
  // Effect of crosslinker on pore size; BMWorld::initializeCaAlg
  BMWorld::PoreSize[0] = static_cast<float>(bm.pore_size.crosslinker_effect);
  // Baseline pore size; BMWorld::initializeCaAlg
  BMWorld::PoreSize[1] = static_cast<float>(bm.pore_size.baseline);
  // Baseline hydrogel mass loss; BMWorld::initializeCaAlg
  BMWorld::MassLoss[0] = static_cast<float>(bm.mass_loss.baseline);
  // Effect of crosslinker on mass loss; BMWorld::initializeCaAlg
  BMWorld::MassLoss[1] = static_cast<float>(bm.mass_loss.crosslinker_effect);
  // Effect of time on mass loss; BMWorld::initializeCaAlg
  BMWorld::MassLoss[2] = static_cast<float>(bm.mass_loss.time_effect);
  // Crosslinker × time interaction on mass loss; BMWorld::initializeCaAlg
  BMWorld::MassLoss[3] =
      static_cast<float>(bm.mass_loss.crosslinker_time_interaction);
  // Baseline swelling ratio; BMWorld::initializeCaAlg
  BMWorld::SwellRatio[0] = static_cast<float>(bm.swell_ratio.baseline);
  // Effect of time on swelling ratio; BMWorld::initializeCaAlg
  BMWorld::SwellRatio[1] = static_cast<float>(bm.swell_ratio.time_effect);
  // Effect of alginate concentration on swelling; BMWorld::initializeCaAlg
  BMWorld::SwellRatio[2] =
      static_cast<float>(bm.swell_ratio.alginate_concentration_effect);
  // Time × crosslinker interaction on swelling; BMWorld::initializeCaAlg
  BMWorld::SwellRatio[3] =
      static_cast<float>(bm.swell_ratio.time_crosslinker_interaction);
  // Alginate × crosslinker interaction on swelling; BMWorld::initializeCaAlg
  BMWorld::SwellRatio[4] =
      static_cast<float>(bm.swell_ratio.alginate_crosslinker_interaction);
}

void log_biology_parameters(const BiologyParametersConfig &cfg,
                            const std::string &source_path)
{
  std::cout << "[biology_parameters] source: " << source_path << std::endl;

  const auto &cp = cfg.cell.proliferation;
  log_entry("cell.proliferation.hours_between_proliferation",
            cp.hours_between_proliferation);
  log_entry("cell.proliferation.tgf_threshold", cp.tgf_threshold);
  log_entry("cell.proliferation.log_scale", cp.log_scale);
  log_entry("cell.proliferation.log_offset", cp.log_offset);

  const auto &cs = cfg.cell.cytokine_synthesis;
  log_entry("cell.cytokine_synthesis.tgf_baseline", cs.tgf_baseline);
  log_entry("cell.cytokine_synthesis.tgf_feedback_tgf", cs.tgf_feedback_tgf);
  log_entry("cell.cytokine_synthesis.tgf_feedback_il1beta",
            cs.tgf_feedback_il1beta);
  log_entry("cell.cytokine_synthesis.tgf_feedback_tnf", cs.tgf_feedback_tnf);
  log_entry("cell.cytokine_synthesis.tnf_baseline", cs.tnf_baseline);
  log_entry("cell.cytokine_synthesis.tnf_feedback_il1beta",
            cs.tnf_feedback_il1beta);
  log_entry("cell.cytokine_synthesis.tnf_feedback_tgf_denom",
            cs.tnf_feedback_tgf_denom);
  log_entry("cell.cytokine_synthesis.il1beta_baseline", cs.il1beta_baseline);
  log_entry("cell.cytokine_synthesis.il1beta_feedback_tnf",
            cs.il1beta_feedback_tnf);
  log_entry("cell.cytokine_synthesis.il1beta_feedback_tgf_denom",
            cs.il1beta_feedback_tgf_denom);

  const auto &st = cfg.stem;
  log_entry("stem.ocr_fmol_per_hour_per_cell", st.ocr_fmol_per_hour_per_cell);
  log_entry("stem.apoptosis_chance", st.apoptosis_chance);
  log_entry("stem.migration.elasticity_effect", st.migration.elasticity_effect);
  log_entry("stem.migration.baseline_speed", st.migration.baseline_speed);
  log_entry("stem.cytokine_synthesis.tgf_baseline",
            st.cytokine_synthesis.tgf_baseline);
  log_entry("stem.cytokine_synthesis.tnf_baseline",
            st.cytokine_synthesis.tnf_baseline);
  log_entry("stem.cytokine_synthesis.il1beta_baseline",
            st.cytokine_synthesis.il1beta_baseline);
  log_entry("stem.collagen_synthesis.baseline_rate",
            st.collagen_synthesis_baseline_rate);
  log_entry("stem.aggrecan_synthesis.tgf_threshold",
            st.aggrecan_synthesis_tgf_threshold);
  log_entry("stem.proliferation.tgf_threshold", st.proliferation.tgf_threshold);
  log_entry("stem.proliferation.tnf_effect", st.proliferation.tnf_effect);
  log_entry("stem.proliferation.il1beta_effect", st.proliferation.il1beta_effect);
  log_entry("stem.proliferation.elasticity_effect",
            st.proliferation.elasticity_effect);
  log_entry("stem.differentiation.asymmetric_probability",
            st.differentiation.asymmetric_probability);
  log_entry("stem.differentiation.baseline_probability",
            st.differentiation.baseline_probability);
  log_entry("stem.differentiation.tgf_effect", st.differentiation.tgf_effect);
  log_entry("stem.differentiation.hours_between_attempts",
            st.differentiation.hours_between_attempts);

  const auto &pg = cfg.progen;
  log_entry("progen.ocr_fmol_per_hour_per_cell", pg.ocr_fmol_per_hour_per_cell);
  log_entry("progen.apoptosis_chance", pg.apoptosis_chance);
  log_entry("progen.migration.elasticity_effect", pg.migration.elasticity_effect);
  log_entry("progen.migration.baseline_speed", pg.migration.baseline_speed);
  log_entry("progen.cytokine_synthesis.tgf_baseline",
            pg.cytokine_synthesis.tgf_baseline);
  log_entry("progen.cytokine_synthesis.tnf_baseline",
            pg.cytokine_synthesis.tnf_baseline);
  log_entry("progen.cytokine_synthesis.il1beta_baseline",
            pg.cytokine_synthesis.il1beta_baseline);
  log_entry("progen.aggrecan_synthesis.baseline_rate",
            pg.aggrecan_synthesis_baseline_rate);

  const auto &np = cfg.np;
  log_entry("np.ocr_fmol_per_hour_per_cell", np.ocr_fmol_per_hour_per_cell);
  log_entry("np.apoptosis_chance", np.apoptosis_chance);
  log_entry("np.migration.elasticity_effect", np.migration.elasticity_effect);
  log_entry("np.migration.baseline_speed", np.migration.baseline_speed);
  log_entry("np.collagen_synthesis.scaling_factor",
            np.collagen_synthesis.scaling_factor);
  log_entry("np.collagen_synthesis.time_effect", np.collagen_synthesis.time_effect);
  log_entry("np.collagen_synthesis.baseline_rate",
            np.collagen_synthesis.baseline_rate);
  log_entry("np.aggrecan_synthesis.scaling_factor",
            np.aggrecan_synthesis.scaling_factor);
  log_entry("np.aggrecan_synthesis.time_effect", np.aggrecan_synthesis.time_effect);
  log_entry("np.aggrecan_synthesis.baseline_rate",
            np.aggrecan_synthesis.baseline_rate);

  const auto &bm = cfg.biomaterial;
  log_entry("biomaterial.elastic_modulus.intercept", bm.elastic_modulus.intercept);
  log_entry("biomaterial.elastic_modulus.alginate_concentration",
            bm.elastic_modulus.alginate_concentration);
  log_entry("biomaterial.elastic_modulus.crosslinker_density",
            bm.elastic_modulus.crosslinker_density);
  log_entry("biomaterial.elastic_modulus.alginate_molecular_weight",
            bm.elastic_modulus.alginate_molecular_weight);
  log_entry("biomaterial.elastic_modulus.alginate_crosslinker_interaction",
            bm.elastic_modulus.alginate_crosslinker_interaction);
  log_entry("biomaterial.elastic_modulus.alginate_mw_interaction",
            bm.elastic_modulus.alginate_mw_interaction);
  log_entry("biomaterial.elastic_modulus.mw_crosslinker_interaction",
            bm.elastic_modulus.mw_crosslinker_interaction);
  log_entry("biomaterial.pore_size.crosslinker_effect",
            bm.pore_size.crosslinker_effect);
  log_entry("biomaterial.pore_size.baseline", bm.pore_size.baseline);
  log_entry("biomaterial.mass_loss.baseline", bm.mass_loss.baseline);
  log_entry("biomaterial.mass_loss.crosslinker_effect",
            bm.mass_loss.crosslinker_effect);
  log_entry("biomaterial.mass_loss.time_effect", bm.mass_loss.time_effect);
  log_entry("biomaterial.mass_loss.crosslinker_time_interaction",
            bm.mass_loss.crosslinker_time_interaction);
  log_entry("biomaterial.swell_ratio.baseline", bm.swell_ratio.baseline);
  log_entry("biomaterial.swell_ratio.time_effect", bm.swell_ratio.time_effect);
  log_entry("biomaterial.swell_ratio.alginate_concentration_effect",
            bm.swell_ratio.alginate_concentration_effect);
  log_entry("biomaterial.swell_ratio.time_crosslinker_interaction",
            bm.swell_ratio.time_crosslinker_interaction);
  log_entry("biomaterial.swell_ratio.alginate_crosslinker_interaction",
            bm.swell_ratio.alginate_crosslinker_interaction);
}

void record_biology_parameters_in_run_params(const BiologyParametersConfig &cfg,
                                             const std::string &source_path)
{
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
  auto add = [&](const char *key, double value) { values[key] = value; };

  const auto &cp = cfg.cell.proliferation;
  add("cell.proliferation.hours_between_proliferation",
      cp.hours_between_proliferation);
  add("cell.proliferation.tgf_threshold", cp.tgf_threshold);
  add("cell.proliferation.log_scale", cp.log_scale);
  add("cell.proliferation.log_offset", cp.log_offset);

  const auto &cs = cfg.cell.cytokine_synthesis;
  add("cell.cytokine_synthesis.tgf_baseline", cs.tgf_baseline);
  add("cell.cytokine_synthesis.tgf_feedback_tgf", cs.tgf_feedback_tgf);
  add("cell.cytokine_synthesis.tgf_feedback_il1beta", cs.tgf_feedback_il1beta);
  add("cell.cytokine_synthesis.tgf_feedback_tnf", cs.tgf_feedback_tnf);
  add("cell.cytokine_synthesis.tnf_baseline", cs.tnf_baseline);
  add("cell.cytokine_synthesis.tnf_feedback_il1beta", cs.tnf_feedback_il1beta);
  add("cell.cytokine_synthesis.tnf_feedback_tgf_denom",
      cs.tnf_feedback_tgf_denom);
  add("cell.cytokine_synthesis.il1beta_baseline", cs.il1beta_baseline);
  add("cell.cytokine_synthesis.il1beta_feedback_tnf", cs.il1beta_feedback_tnf);
  add("cell.cytokine_synthesis.il1beta_feedback_tgf_denom",
      cs.il1beta_feedback_tgf_denom);

  const auto &st = cfg.stem;
  add("stem.ocr_fmol_per_hour_per_cell", st.ocr_fmol_per_hour_per_cell);
  add("stem.apoptosis_chance", st.apoptosis_chance);
  add("stem.migration.elasticity_effect", st.migration.elasticity_effect);
  add("stem.migration.baseline_speed", st.migration.baseline_speed);
  add("stem.cytokine_synthesis.tgf_baseline", st.cytokine_synthesis.tgf_baseline);
  add("stem.cytokine_synthesis.tnf_baseline", st.cytokine_synthesis.tnf_baseline);
  add("stem.cytokine_synthesis.il1beta_baseline",
      st.cytokine_synthesis.il1beta_baseline);
  add("stem.collagen_synthesis.baseline_rate", st.collagen_synthesis_baseline_rate);
  add("stem.aggrecan_synthesis.tgf_threshold",
      st.aggrecan_synthesis_tgf_threshold);
  add("stem.proliferation.tgf_threshold", st.proliferation.tgf_threshold);
  add("stem.proliferation.tnf_effect", st.proliferation.tnf_effect);
  add("stem.proliferation.il1beta_effect", st.proliferation.il1beta_effect);
  add("stem.proliferation.elasticity_effect", st.proliferation.elasticity_effect);
  add("stem.differentiation.asymmetric_probability",
      st.differentiation.asymmetric_probability);
  add("stem.differentiation.baseline_probability",
      st.differentiation.baseline_probability);
  add("stem.differentiation.tgf_effect", st.differentiation.tgf_effect);
  add("stem.differentiation.hours_between_attempts",
      st.differentiation.hours_between_attempts);

  const auto &pg = cfg.progen;
  add("progen.ocr_fmol_per_hour_per_cell", pg.ocr_fmol_per_hour_per_cell);
  add("progen.apoptosis_chance", pg.apoptosis_chance);
  add("progen.migration.elasticity_effect", pg.migration.elasticity_effect);
  add("progen.migration.baseline_speed", pg.migration.baseline_speed);
  add("progen.cytokine_synthesis.tgf_baseline", pg.cytokine_synthesis.tgf_baseline);
  add("progen.cytokine_synthesis.tnf_baseline", pg.cytokine_synthesis.tnf_baseline);
  add("progen.cytokine_synthesis.il1beta_baseline",
      pg.cytokine_synthesis.il1beta_baseline);
  add("progen.aggrecan_synthesis.baseline_rate",
      pg.aggrecan_synthesis_baseline_rate);

  const auto &np = cfg.np;
  add("np.ocr_fmol_per_hour_per_cell", np.ocr_fmol_per_hour_per_cell);
  add("np.apoptosis_chance", np.apoptosis_chance);
  add("np.migration.elasticity_effect", np.migration.elasticity_effect);
  add("np.migration.baseline_speed", np.migration.baseline_speed);
  add("np.collagen_synthesis.scaling_factor", np.collagen_synthesis.scaling_factor);
  add("np.collagen_synthesis.time_effect", np.collagen_synthesis.time_effect);
  add("np.collagen_synthesis.baseline_rate", np.collagen_synthesis.baseline_rate);
  add("np.aggrecan_synthesis.scaling_factor", np.aggrecan_synthesis.scaling_factor);
  add("np.aggrecan_synthesis.time_effect", np.aggrecan_synthesis.time_effect);
  add("np.aggrecan_synthesis.baseline_rate", np.aggrecan_synthesis.baseline_rate);

  const auto &bm = cfg.biomaterial;
  add("biomaterial.elastic_modulus.intercept", bm.elastic_modulus.intercept);
  add("biomaterial.elastic_modulus.alginate_concentration",
      bm.elastic_modulus.alginate_concentration);
  add("biomaterial.elastic_modulus.crosslinker_density",
      bm.elastic_modulus.crosslinker_density);
  add("biomaterial.elastic_modulus.alginate_molecular_weight",
      bm.elastic_modulus.alginate_molecular_weight);
  add("biomaterial.elastic_modulus.alginate_crosslinker_interaction",
      bm.elastic_modulus.alginate_crosslinker_interaction);
  add("biomaterial.elastic_modulus.alginate_mw_interaction",
      bm.elastic_modulus.alginate_mw_interaction);
  add("biomaterial.elastic_modulus.mw_crosslinker_interaction",
      bm.elastic_modulus.mw_crosslinker_interaction);
  add("biomaterial.pore_size.crosslinker_effect", bm.pore_size.crosslinker_effect);
  add("biomaterial.pore_size.baseline", bm.pore_size.baseline);
  add("biomaterial.mass_loss.baseline", bm.mass_loss.baseline);
  add("biomaterial.mass_loss.crosslinker_effect", bm.mass_loss.crosslinker_effect);
  add("biomaterial.mass_loss.time_effect", bm.mass_loss.time_effect);
  add("biomaterial.mass_loss.crosslinker_time_interaction",
      bm.mass_loss.crosslinker_time_interaction);
  add("biomaterial.swell_ratio.baseline", bm.swell_ratio.baseline);
  add("biomaterial.swell_ratio.time_effect", bm.swell_ratio.time_effect);
  add("biomaterial.swell_ratio.alginate_concentration_effect",
      bm.swell_ratio.alginate_concentration_effect);
  add("biomaterial.swell_ratio.time_crosslinker_interaction",
      bm.swell_ratio.time_crosslinker_interaction);
  add("biomaterial.swell_ratio.alginate_crosslinker_interaction",
      bm.swell_ratio.alginate_crosslinker_interaction);

  if (!root.contains("simulation"))
    root["simulation"] = json::object();
  root["simulation"]["simulation_config"] = source_path;
  root["biology_parameters"] = values;

  std::ofstream out(path);
  if (!out)
  {
    std::fprintf(stderr, "Warning: cannot update run params at %s\n", path);
    return;
  }
  out << root.dump(2) << std::endl;
}
