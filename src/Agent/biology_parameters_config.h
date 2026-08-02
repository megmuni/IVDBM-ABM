#ifndef IVDBM_BIOLOGY_PARAMETERS_CONFIG_H
#define IVDBM_BIOLOGY_PARAMETERS_CONFIG_H

/**
 * @file biology_parameters_config.h
 * @brief JSON-backed cell rule and hydrogel calibration parameters.
 *
 * Loaded from the @c biology section of simulation_config.json.
 */

#include <map>
#include <string>
#include <vector>

struct CellProliferationParams {
  double hours_between_proliferation = 0;
  double tgf_threshold = 0;
  double log_scale = 0;
  double log_offset = 0;
};

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

struct MigrationParams {
  double elasticity_effect = 0;
  double baseline_speed = 0;
};

struct StemCytokineSynthesisParams {
  double tgf_baseline = 0;
  double tnf_baseline = 0;
  double il1beta_baseline = 0;
};

struct ProgenCytokineSynthesisParams {
  double tgf_baseline = 0;
  double tnf_baseline = 0;
  double il1beta_baseline = 0;
};

struct StemProliferationParams {
  double tgf_threshold = 0;
  double tnf_effect = 0;
  double il1beta_effect = 0;
  double elasticity_effect = 0;
};

struct StemDifferentiationParams {
  double asymmetric_probability = 0;
  double baseline_probability = 0;
  double tgf_effect = 0;
  /** Hours between differentiation attempts. */
  double hours_between_attempts = 0;
};

struct NpCollagenSynthesisParams {
  double scaling_factor = 0;
  double time_effect = 0;
  double baseline_rate = 0;
};

struct NpAggrecanSynthesisParams {
  double scaling_factor = 0;
  double time_effect = 0;
  double baseline_rate = 0;
};

struct ElasticModulusParams {
  double intercept = 0;
  double alginate_concentration = 0;
  double crosslinker_density = 0;
  double alginate_molecular_weight = 0;
  double alginate_crosslinker_interaction = 0;
  double alginate_mw_interaction = 0;
  double mw_crosslinker_interaction = 0;
};

struct PoreSizeParams {
  double crosslinker_effect = 0;
  double baseline = 0;
};

struct MassLossParams {
  double baseline = 0;
  double crosslinker_effect = 0;
  double time_effect = 0;
  double crosslinker_time_interaction = 0;
};

struct SwellRatioParams {
  double baseline = 0;
  double time_effect = 0;
  double alginate_concentration_effect = 0;
  double time_crosslinker_interaction = 0;
  double alginate_crosslinker_interaction = 0;
};

struct CellParams {
  CellProliferationParams proliferation;
  CellCytokineSynthesisParams cytokine_synthesis;
};

struct StemParams {
  /** OCR in fmol/h/cell (converted to per-tick in apply). */
  double ocr_fmol_per_hour_per_cell = 0;
  double apoptosis_chance = 0;
  MigrationParams migration;
  StemCytokineSynthesisParams cytokine_synthesis;
  double collagen_synthesis_baseline_rate = 0;
  double aggrecan_synthesis_tgf_threshold = 0;
  StemProliferationParams proliferation;
  StemDifferentiationParams differentiation;
};

struct ProgenParams {
  double ocr_fmol_per_hour_per_cell = 0;
  double apoptosis_chance = 0;
  MigrationParams migration;
  ProgenCytokineSynthesisParams cytokine_synthesis;
  double aggrecan_synthesis_baseline_rate = 0;
};

struct NpParams {
  double ocr_fmol_per_hour_per_cell = 0;
  double apoptosis_chance = 0;
  MigrationParams migration;
  NpCollagenSynthesisParams collagen_synthesis;
  NpAggrecanSynthesisParams aggrecan_synthesis;
};

struct BiomaterialParams {
  ElasticModulusParams elastic_modulus;
  PoreSizeParams pore_size;
  MassLossParams mass_loss;
  SwellRatioParams swell_ratio;
};

struct BiologyParametersConfig {
  CellParams cell;
  StemParams stem;
  ProgenParams progen;
  NpParams np;
  BiomaterialParams biomaterial;
};

BiologyParametersConfig load_biology_parameters_config(const std::string &path);

void apply_biology_parameters(const BiologyParametersConfig &cfg);

void log_biology_parameters(const BiologyParametersConfig &cfg,
                            const std::string &source_path);

void record_biology_parameters_in_run_params(const BiologyParametersConfig &cfg,
                                             const std::string &source_path);

#endif
