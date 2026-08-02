#ifndef IVDBM_BIOLOGY_PARAMETERS_CONFIG_H
#define IVDBM_BIOLOGY_PARAMETERS_CONFIG_H

/**
 * @file biology_parameters_config.h
 * @brief JSON-backed cell rule and hydrogel calibration parameters.
 *
 * Loaded from the @c biology section of simulation_config.json.
 */

#include <string>

#include <nlohmann/json.hpp>

struct CellProliferationParams {
  double hours_between_proliferation = 0;
  double tgf_threshold = 0;
  double log_scale = 0;
  double log_offset = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CellProliferationParams,
                                   hours_between_proliferation, tgf_threshold,
                                   log_scale, log_offset)

struct CellCytokineSynthesisParams {
  double tgf_baseline = 0;
  double tgf_feedback_tgf = 0;
  double tgf_feedback_il1beta = 0;
  double tgf_feedback_tnf = 0;
  double tnf_baseline = 0;
  double tnf_feedback_il1beta = 0;
  double tnf_feedback_tgf_denom = 0;
  double il1beta_baseline = 0;
  double il1beta_feedback_tnf = 0;
  double il1beta_feedback_tgf_denom = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CellCytokineSynthesisParams, tgf_baseline,
                                   tgf_feedback_tgf, tgf_feedback_il1beta,
                                   tgf_feedback_tnf, tnf_baseline,
                                   tnf_feedback_il1beta, tnf_feedback_tgf_denom,
                                   il1beta_baseline, il1beta_feedback_tnf,
                                   il1beta_feedback_tgf_denom)

struct MigrationParams {
  double elasticity_effect = 0;
  double baseline_speed = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MigrationParams, elasticity_effect,
                                   baseline_speed)

struct StemCytokineSynthesisParams {
  double tgf_baseline = 0;
  double tnf_baseline = 0;
  double il1beta_baseline = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StemCytokineSynthesisParams, tgf_baseline,
                                   tnf_baseline, il1beta_baseline)

struct ProgenCytokineSynthesisParams {
  double tgf_baseline = 0;
  double tnf_baseline = 0;
  double il1beta_baseline = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ProgenCytokineSynthesisParams, tgf_baseline,
                                   tnf_baseline, il1beta_baseline)

struct StemProliferationParams {
  double tgf_threshold = 0;
  double tnf_effect = 0;
  double il1beta_effect = 0;
  double elasticity_effect = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StemProliferationParams, tgf_threshold,
                                   tnf_effect, il1beta_effect,
                                   elasticity_effect)

struct StemDifferentiationParams {
  double asymmetric_probability = 0;
  double baseline_probability = 0;
  double tgf_effect = 0;
  /** Hours between differentiation attempts. */
  double hours_between_attempts = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StemDifferentiationParams,
                                   asymmetric_probability, baseline_probability,
                                   tgf_effect, hours_between_attempts)

struct StemCollagenSynthesisParams {
  double baseline_rate = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StemCollagenSynthesisParams, baseline_rate)

struct StemAggrecanSynthesisParams {
  double tgf_threshold = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StemAggrecanSynthesisParams, tgf_threshold)

struct ProgenAggrecanSynthesisParams {
  double baseline_rate = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ProgenAggrecanSynthesisParams, baseline_rate)

struct NpCollagenSynthesisParams {
  double scaling_factor = 0;
  double time_effect = 0;
  double baseline_rate = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NpCollagenSynthesisParams, scaling_factor,
                                   time_effect, baseline_rate)

struct NpAggrecanSynthesisParams {
  double scaling_factor = 0;
  double time_effect = 0;
  double baseline_rate = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NpAggrecanSynthesisParams, scaling_factor,
                                   time_effect, baseline_rate)

struct ElasticModulusParams {
  double intercept = 0;
  double alginate_concentration = 0;
  double crosslinker_density = 0;
  double alginate_molecular_weight = 0;
  double alginate_crosslinker_interaction = 0;
  double alginate_mw_interaction = 0;
  double mw_crosslinker_interaction = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ElasticModulusParams, intercept,
                                   alginate_concentration, crosslinker_density,
                                   alginate_molecular_weight,
                                   alginate_crosslinker_interaction,
                                   alginate_mw_interaction,
                                   mw_crosslinker_interaction)

struct PoreSizeParams {
  double crosslinker_effect = 0;
  double baseline = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(PoreSizeParams, crosslinker_effect, baseline)

struct MassLossParams {
  double baseline = 0;
  double crosslinker_effect = 0;
  double time_effect = 0;
  double crosslinker_time_interaction = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MassLossParams, baseline, crosslinker_effect,
                                   time_effect, crosslinker_time_interaction)

struct SwellRatioParams {
  double baseline = 0;
  double time_effect = 0;
  double alginate_concentration_effect = 0;
  double time_crosslinker_interaction = 0;
  double alginate_crosslinker_interaction = 0;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SwellRatioParams, baseline, time_effect,
                                   alginate_concentration_effect,
                                   time_crosslinker_interaction,
                                   alginate_crosslinker_interaction)

struct CellParams {
  CellProliferationParams proliferation;
  CellCytokineSynthesisParams cytokine_synthesis;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(CellParams, proliferation,
                                   cytokine_synthesis)

struct StemParams {
  /** OCR in fmol/h/cell (converted to per-tick in apply). */
  double ocr_fmol_per_hour_per_cell = 0;
  double apoptosis_chance = 0;
  MigrationParams migration;
  StemCytokineSynthesisParams cytokine_synthesis;
  StemCollagenSynthesisParams collagen_synthesis;
  StemAggrecanSynthesisParams aggrecan_synthesis;
  StemProliferationParams proliferation;
  StemDifferentiationParams differentiation;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(StemParams, ocr_fmol_per_hour_per_cell,
                                   apoptosis_chance, migration,
                                   cytokine_synthesis, collagen_synthesis,
                                   aggrecan_synthesis, proliferation,
                                   differentiation)

struct ProgenParams {
  double ocr_fmol_per_hour_per_cell = 0;
  double apoptosis_chance = 0;
  MigrationParams migration;
  ProgenCytokineSynthesisParams cytokine_synthesis;
  ProgenAggrecanSynthesisParams aggrecan_synthesis;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ProgenParams, ocr_fmol_per_hour_per_cell,
                                   apoptosis_chance, migration,
                                   cytokine_synthesis, aggrecan_synthesis)

struct NpParams {
  double ocr_fmol_per_hour_per_cell = 0;
  double apoptosis_chance = 0;
  MigrationParams migration;
  NpCollagenSynthesisParams collagen_synthesis;
  NpAggrecanSynthesisParams aggrecan_synthesis;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NpParams, ocr_fmol_per_hour_per_cell,
                                   apoptosis_chance, migration,
                                   collagen_synthesis, aggrecan_synthesis)

struct BiomaterialParams {
  ElasticModulusParams elastic_modulus;
  PoreSizeParams pore_size;
  MassLossParams mass_loss;
  SwellRatioParams swell_ratio;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BiomaterialParams, elastic_modulus,
                                   pore_size, mass_loss, swell_ratio)

struct BiologyParametersConfig {
  CellParams cell;
  StemParams stem;
  ProgenParams progen;
  NpParams np;
  BiomaterialParams biomaterial;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BiologyParametersConfig, cell, stem, progen,
                                   np, biomaterial)

BiologyParametersConfig load_biology_parameters_config(const std::string &path);

void apply_biology_parameters(const BiologyParametersConfig &cfg);

void log_biology_parameters(const BiologyParametersConfig &cfg,
                            const std::string &source_path);

void record_biology_parameters_in_run_params(const BiologyParametersConfig &cfg,
                                             const std::string &source_path);

#endif
