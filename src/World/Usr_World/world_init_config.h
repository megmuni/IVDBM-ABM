#ifndef IVDBM_WORLD_INIT_CONFIG_H
#define IVDBM_WORLD_INIT_CONFIG_H

/**
 * @file world_init_config.h
 * @brief JSON-backed scaffold seeding / alginate hydrogel composition
 * parameters, loaded from the @c world_init section of simulation_config.json
 *
 */

#include <string>

struct WorldInitAlginateParams {
  double wv_percent = 0.0;    // Alginate concentration (% w/v)
  double high_mw_ratio = 0.0; // Ratio component of high MW alginate
  double low_mw_ratio = 0.0;  // Ratio component of low MW alginate
  double ca_mm = 0.0;         // Calcium crosslinker concentration (mM)
};

struct WorldInitParams {
  /**
   * Initial MSC seed count for the scaffold; 0 selects the default seeding
   * density (see BMWorld::initializeCells). Consumed by BMWorld::userInput.
   */
  int msc_count = 0;
  WorldInitAlginateParams alginate;
  /** Peptide conjugation type; only used when compiled with PEPTIDE_BM. */
  std::string peptide;
};

/**
 * @brief Load and validate the @c world_init section of a simulation config.
 * @param path Path to the JSON config file (simulation_config.json).
 */
WorldInitParams load_world_init_config(const std::string &path);

#endif
