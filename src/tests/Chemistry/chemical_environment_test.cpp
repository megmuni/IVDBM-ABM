#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../../common.h"
#include "../../enums.h"
#include "chemical_environment.h"
#include "species_registry.h"

#ifndef IVDBM_CONFIG_DIR
#define IVDBM_CONFIG_DIR "configFiles"
#endif

namespace {

std::string test_simulation_config_path()
{
    return std::string(IVDBM_CONFIG_DIR) + "/simulation_config.template.json";
}

void load_test_env(ChemicalEnvironment &env)
{
    env.load_from_config(test_simulation_config_path(), 1.0, 1.0);
    env.allocate_channels_from_config();
}

} // namespace

TEST_CASE("ChemicalEnvironment merge updates concentration from delta channel",
          "[chemistry][environment]") {
  ChemicalEnvironment env(2, 2, 2);
  load_test_env(env);

  env.set_concentration(0, TNF, 10.f);
  env.accumulate_secretion(0, TNF, 2.f);

  env.merge_and_reset_secretion();

  REQUIRE(env.concentration_at(0, TNF) == Approx(12.f));
  REQUIRE(env.channels(TNF).secretion_delta[0] == Approx(0.f));
}

TEST_CASE("ChemicalEnvironment merge clamps concentration to non-negative",
          "[chemistry][environment]") {
  ChemicalEnvironment env(2, 2, 2);
  load_test_env(env);

  env.set_concentration(0, TNF, 5.f);
  env.accumulate_secretion(0, TNF, -12.f);
  env.merge_and_reset_secretion();

  REQUIRE(env.concentration_at(0, TNF) == Approx(0.f));
}

TEST_CASE("ChemicalEnvironment agent API reads and accumulates secretion",
          "[chemistry][environment][agent]") {
  ChemicalEnvironment env(2, 2, 2);
  load_test_env(env);

  env.set_concentration(3, TGF, 7.f);
  REQUIRE(env.concentration_at(3, TGF) == Approx(7.f));
  REQUIRE(env.concentration_at_channel(3, pTGF) == Approx(7.f));

  env.accumulate_secretion(3, TNF, 1.5f);
  REQUIRE(env.channels(TNF).secretion_delta[3] == Approx(1.5f));

  env.set_concentration(2, TGF, 3.5f);
  env.update_chemotaxis_from_species(TGF);
  REQUIRE(env.chemotaxis_at(2) == Approx(3.5f));
}

TEST_CASE("ChemicalEnvironment baseline helpers set concentration and totals",
          "[chemistry][environment][baseline]") {
  ChemicalEnvironment env(2, 2, 1);
  load_test_env(env);

  env.set_concentration(1, TGF, 5.f);
  REQUIRE(env.concentration_at(1, TGF) == Approx(5.f));

  env.update_chemotaxis_from_species(TGF);
  REQUIRE(env.chemotaxis_at(1) == Approx(5.f));

  env.recompute_world_totals();
  REQUIRE(env.total_tgf() == Approx(5.f));
}

TEST_CASE("ChemicalEnvironment loads tick interval from JSON",
          "[chemistry][environment][config]") {
  ChemicalEnvironment env(2, 2, 2);
  load_test_env(env);
  REQUIRE(env.tick_interval_minutes() == Approx(30.0));
  REQUIRE(env.baseline_total_mass_for("TNF") == Approx(0.f));
}
